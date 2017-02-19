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
  if ( monitor1 != NULL ) delete monitor1;
  if ( monitor2 != NULL ) delete monitor2;
}

void CoccolithSimulation::loadVoxels( const char* fname )
{
  material.loadRaw( fname );
  cout << material.sizeX() << " " << material.sizeY() << " " << material.sizeZ() << endl;
  cout << getPMLThickness() << endl;
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
      crn1 = meep::vec( getSrcPos(), 0, 0 );
      crn2 = meep::vec( getSrcPos(), material.sizeY(), material.sizeZ() );
      break;
    case MainPropDirection_t::Y:
      crn1 = meep::vec( 0, getSrcPos(), 0 );
      crn2 = meep::vec( material.sizeX(), getSrcPos(), material.sizeZ() );
      break;
    case MainPropDirection_t::Z:
      crn1 = meep::vec( 0, 0, getSrcPos() );
      crn2 = meep::vec( material.sizeX(), material.sizeY(), getSrcPos() );
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
  domainInfo();
  addSourceVolume();
  addStructure();
  addFields();
  addSource();
  addFluxPlanes();
  isInitialized = true;
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

  meep::vec sourceDFTCrn1 = crn1 + displacement*getSrcFluxPos();
  meep::vec sourceDFTCrn2 = crn2 + displacement*getSrcFluxPos();

  if ( dftVolSource != NULL ) delete dftVolSource;
  dftVolSource = new meep::volume( sourceDFTCrn1, sourceDFTCrn2 );

  meep::vec domainDiagonal( material.sizeX(), material.sizeY(), material.sizeZ() );
  double size = domainDiagonal&displacement; // In MEEP: & is dot product

  meep::vec transDFTCrn1 = crn1 + displacement*( getTransFluxPos() - pmlThicknessInWavelengths*getWavelength() );
  meep::vec transDFTCrn2 = crn2 + displacement*( getTransFluxPos() - pmlThicknessInWavelengths*getWavelength() );

  if ( dftVolTransmit != NULL ) delete dftVolTransmit;
  dftVolTransmit = new meep::volume( transDFTCrn1, transDFTCrn2 );

  if ( srcFlux != NULL ) delete srcFlux;
  if ( transmitFlux != NULL ) delete transmitFlux;
  srcFlux = new meep::dft_flux( field->add_dft_flux_plane( *dftVolSource, centerFrequency-freqWidth/2.0, centerFrequency+freqWidth/2.0, nfreq ) );
  transmitFlux = new meep::dft_flux( field->add_dft_flux_plane( *dftVolTransmit, centerFrequency-freqWidth/2.0, centerFrequency+freqWidth/2.0, nfreq ) );
}


double CoccolithSimulation::getLowerBorderInPropDir() const
{
  return pmlThicknessInWavelengths*getWavelength();
}

double CoccolithSimulation::getUpperBorderInPropDir() const
{
  switch( propagationDir )
  {
    case MainPropDirection_t::X:
      return material.sizeX() - pmlThicknessInWavelengths*getWavelength();
    case MainPropDirection_t::Y:
      return material.sizeY() - pmlThicknessInWavelengths*getWavelength();
    case MainPropDirection_t::Z:
      return material.sizeZ() - pmlThicknessInWavelengths*getWavelength();
  }
}
double CoccolithSimulation::getSrcPos() const
{
  switch( srcPos )
  {
    case SourcePosition_t::TOP:
      return getLowerBorderInPropDir();
    case SourcePosition_t::BOTTOM:
      return getUpperBorderInPropDir();
  }
}

double CoccolithSimulation::getSrcFluxPos() const
{
  switch( srcPos )
  {
    case SourcePosition_t::TOP:
      return getSrcPos()+5.0;
    case SourcePosition_t::BOTTOM:
      return getSrcPos() - 5.0;
  }
}

double CoccolithSimulation::getTransFluxPos() const
{
  switch( srcPos )
  {
    case SourcePosition_t::TOP:
      return getUpperBorderInPropDir();
    case SourcePosition_t::BOTTOM:
      return getLowerBorderInPropDir();
  }
}

void CoccolithSimulation::setMonitorPlanes()
{
  if ( monitor1 != NULL ) delete monitor1;
  if ( monitor2 != NULL ) delete monitor2;

  meep::vec vMain;
  meep::vec v1;
  meep::vec v2;
  meep::vec orig1;
  meep::vec orig2;
  double dx = static_cast<double>(material.sizeX())/static_cast<double>(nMonitorX);
  double dy = static_cast<double>(material.sizeY())/static_cast<double>(nMonitorY);
  double dz = static_cast<double>(material.sizeZ())/static_cast<double>(nMonitorZ);
  switch ( propagationDir )
  {
    case MainPropDirection_t::X:
      vMain = meep::vec(dx, 0.0, 0.0 );
      v1 = meep::vec(0.0, dy, 0.0 );
      v2 = meep::vec(0.0, 0.0, dz );
      monitor1 = new FieldMonitor( nMonitorY, nMonitorX );
      monitor2 = new FieldMonitor( nMonitorZ, nMonitorX );
      monitor1->setName( "XY-plane" );
      monitor2->setName( "XZ-plane" );
      orig1 = meep::vec( 0.0, 0.0, material.sizeZ()/2.0 );
      orig2 = meep::vec( 0.0, material.sizeY()/2.0, 0.0 );
      break;
    case MainPropDirection_t::Y:
      vMain = meep::vec(0.0, dy, 0.0 );
      v1 = meep::vec( dx, 0.0, 0.0 );
      v2 = meep::vec( 0.0, 0.0, dz );
      monitor1 = new FieldMonitor( nMonitorX, nMonitorY );
      monitor2 = new FieldMonitor( nMonitorZ, nMonitorY );
      monitor1->setName( "XY-plane" );
      monitor2->setName( "YZ-plane" );
      orig1 = meep::vec( 0.0, 0.0, material.sizeZ()/2.0 );
      orig2 = meep::vec( material.sizeX()/2.0, 0.0, 0.0 );
      break;
    case MainPropDirection_t::Z:
      vMain = meep::vec(0.0, 0.0, dz );
      v1 = meep::vec( dx, 0.0, 0.0 );
      v2 = meep::vec( 0.0, dy, 0.0 );
      monitor1 = new FieldMonitor( nMonitorX, nMonitorZ );
      monitor2 = new FieldMonitor( nMonitorY, nMonitorZ );
      monitor1->setName( "XZ-plane" );
      monitor2->setName( "YZ-plane" );
      orig1 = meep::vec( 0.0, material.sizeY()/2.0, 0.0 );
      orig2 = meep::vec( material.sizeX()/2.0, 0.0, 0.0 );
      break;
  }

  monitor1->setDisplacementVectors( v1, vMain );
  monitor2->setDisplacementVectors( v2, vMain );
  monitor1->setOrigin( orig1 );
  monitor2->setOrigin( orig2 );

  // Add plots
  plots.addPlot( monitor1->getName().c_str() );
  plots.addPlot( monitor2->getName().c_str() );
}

void CoccolithSimulation::run()
{
  if ( !isInitialized )
  {
    throw( runtime_error("Simulation is not initialized!") );
  }

  unsigned int iter = 0;
  while ( field->time() < tEnd )
  {
    field->step();

    if ( iter%plotUpdateFreq == 0 )
    {
      visualize();
    }
    iter++;
  }
}

void CoccolithSimulation::visualize()
{
  monitor1->setIntensity( *field );
  monitor2->setIntensity( *field );

  plots.get( monitor1->getName().c_str() ).fillVertexArray( monitor1->get() );
  plots.get( monitor2->getName().c_str() ).fillVertexArray( monitor2->get() );
  plots.show();
}

void CoccolithSimulation::domainInfo() const
{
  double speedOfLight = 2.997E8; // nm/ns
  double L = material.getVoxelSize();
  double PI = acos(-1.0);
  cout << "-----------------------------------------------------------\n";
  cout << "Domain info:\n";
  cout << "SizeX: " << material.sizeX()*L << " nm, SizeY: " << material.sizeY()*L << " nm. SizeZ: " << material.sizeZ()*L << " nm\n";
  cout << "dx=dy=dz = " << L/resolution << " nm\n";
  cout << "PML thickness: " << getPMLThickness()*L << " nm\n";
  cout << "Main frequency: " << centerFrequency*speedOfLight/(2.0*PI*L*1000.0) << " (THz)\n";
  cout << "Frequency width: " << freqWidth*speedOfLight/(2.0*PI*L*1000.0) << " (THz)\n";
  cout << "Wavelength: " << getWavelength()*L << " nm\n";
}
