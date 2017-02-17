#include "coccolithSim.hpp"
#include <sstream>
#include <stdexcept>
#include <cmath>

using namespace std;

meep::vec CoccolithSimulation::waveVec;

CoccolithSimulation::~CoccolithSimulation()
{
  if ( srcVol != NULL ) delete srcVol;
  if ( struc != NULL ) delete struc;
  if ( field != NULL ) delete field;
  if ( source != NULL ) delete source;
  if ( dftVolSource != NULL ) delete dftVolSource;
  if ( dftVolTransmit != NULL ) delete dftVolTransmit;
  if ( srcFlux != NULL ) delete srcFlux;
  if ( transmitFlux != NULL ) delete transmitFlux;
}

void CoccolithSimulation::loadVoxels( const char* fname )
{
  material.loadRaw( fname );
  gdvol = meep::vol3d( material.sizeX(), material.sizeY(), material.sizeZ(), resolution );
  materialLoaded = true;
}

cdouble CoccolithSimulation::amplitude( const meep::vec &r )
{
  cdouble im(0.0,1.0);
  return exp( im*(waveVec&r) ); // In MEEP: & is dot product
}

void CoccolithSimulation::setMainPropagationDirection( MainPropDirection_t propDir )
{
  propagationDir = propDir;
  addSourceVolume();
}

void CoccolithSimulation::addSourceVolume()
{
  double pml = getPMLThickness();
  switch( propagationDir )
  {
    case MainPropDirection_t::X:
      crn1 = meep::vec( pml, 0, 0 );
      crn2 = meep::vec( pml, material.sizeY(), material.sizeZ() );
      break;
    case MainPropDirection_t::Y:
      crn1 = meep::vec( 0, pml, 0 );
      crn2 = meep::vec( material.sizeX(), pml, material.sizeZ() );
      break;
    case MainPropDirection_t::Z:
      crn1 = meep::vec( 0, 0, pml );
      crn2 = meep::vec( material.sizeX(), material.sizeY(), pml );
      break;
  }

  if ( srcVol != NULL ) delete srcVol;
  srcVol = new meep::volume( crn1, crn2 );
}

void CoccolithSimulation::addStructure()
{

  if ( source == NULL )
  {
    throw( runtime_error("The source profile needs to be initialized before the structure is added!\n") );
  }
  if ( struc != NULL ) delete struc;

  struc = new meep::structure( gdvol, material, meep::pml( getPMLThickness() ) );
}

double CoccolithSimulation::getWavelength() const
{
  double PI = acos(-1.0);
  return 2*PI/centerFrequency; // In MEEP units c=1
}

void CoccolithSimulation::addFields()
{
  if ( struc == NULL )
  {
    throw( runtime_error("A structure must be added before fields!") );
  }

  if ( field != NULL ) delete field;
  field = new meep::fields( struc );
}

void CoccolithSimulation::initSource( double freq, double fwidth )
{
  centerFrequency = freq;
  freqWidth = fwidth;

  if ( source != NULL ) delete source;
  source = new meep::gaussian_src_time( freq, fwidth );
}

void CoccolithSimulation::addSource()
{
  if ( field == NULL )
  {
    throw( runtime_error("Field needs to be added before a source is added!\n") );
  }
  if ( source == NULL )
  {
    throw( runtime_error("A source profile needs to be added before a source is added! Call initSource!\n") );
  }
  if ( srcVol == NULL )
  {
    throw( runtime_error("A source volume needs to be added before a source is added!\n") );
  }

  meep::component fieldComp;
  switch( propagationDir )
  {
    case MainPropDirection_t::X:
      fieldComp = meep::component::Ey;
      break;
    case MainPropDirection_t::Y:
      fieldComp = meep::component::Ez;
      break;
    case MainPropDirection_t::Z:
      fieldComp = meep::component::Ex;
      break;
  }

  field->add_volume_source( fieldComp, *source, *srcVol, amplitude );
}

void CoccolithSimulation::init()
{
  if ( !materialLoaded )
  {
    throw( runtime_error("No material loaded! Call loadVoxels!") );
  }
  addSourceVolume();
  addStructure();
  addFields();
  addSource();
  addFluxPlanes();
}

double CoccolithSimulation::getPMLThickness() const
{
  return pmlThicknessInWavelengths*getWavelength();
}

void CoccolithSimulation::addFluxPlanes()
{

  if ( source == NULL )
  {
    throw( runtime_error("A source profile needs to be set before flux planes are added!") );
  }
  meep::vec diagonal = crn2 - crn1;
  meep::vec displacement;
  const double ZERO = 1E-12;

  // This vector should be in a plane --> one component is zero
  if (( abs(diagonal.x()) < ZERO ) && ( abs(diagonal.y() > ZERO) ) && ( abs(diagonal.z() > ZERO ) ))
  {
    displacement = meep::vec(1.0, 0, 0 );
  }
  else if (( abs(diagonal.y()) < ZERO ) && ( abs(diagonal.z() > ZERO) ) && ( abs(diagonal.x() > ZERO ) ))
  {
    displacement = meep::vec(0.0, 1.0, 0.0 );
  }
  else if (( abs(diagonal.z()) < ZERO ) && ( abs(diagonal.x() > ZERO) ) && ( abs(diagonal.y() > ZERO ) ))
  {
    displacement = meep::vec( 0.0, 0.0, 1.0  );
  }
  else
  {
    throw ( runtime_error("The diagonal vector obtained from the corners does not lie in a plane orthogonal to one of the coordinate axes!"));
  }

  meep::vec sourceDFTCrn1 = crn1 + displacement*5.0;
  meep::vec sourceDFTCrn2 = crn2 + displacement*5.0;

  if ( dftVolSource != NULL ) delete dftVolSource;
  dftVolSource = new meep::volume( sourceDFTCrn1, sourceDFTCrn2 );

  meep::vec domainDiagonal( material.sizeX(), material.sizeY(), material.sizeZ() );
  double size = domainDiagonal&displacement; // In MEEP: & is dot product

  meep::vec transDFTCrn1 = crn1 + displacement*size - displacement*2.0*getPMLThickness() - displacement*5.0;
  meep::vec transDFTCrn2 = crn2 + displacement*size - displacement*2.0*getPMLThickness() - displacement*5.0;

  if ( dftVolTransmit != NULL ) delete dftVolTransmit;
  dftVolTransmit = new meep::volume( transDFTCrn1, transDFTCrn2 );

  if ( srcFlux != NULL ) delete srcFlux;
  if ( transmitFlux != NULL ) delete transmitFlux;
  srcFlux = new meep::dft_flux( field->add_dft_flux_plane( *dftVolSource, centerFrequency-freqWidth/2.0, centerFrequency+freqWidth/2.0, nfreq ) );
  transmitFlux = new meep::dft_flux( field->add_dft_flux_plane( *dftVolTransmit, centerFrequency-freqWidth/2.0, centerFrequency+freqWidth/2.0, nfreq ) );
}
