#include "coccolithSim.hpp"
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <omp.h>
#include <ctime>
#include <iomanip>
#include <array>
#include <gsl/gsl_integration.h>
#include <mpi.h>
#include <cstring>
#define DEBUG_N2F_INITIALIZATION

using namespace std;

meep::vec CoccolithSimulation::waveVec;
array<Stokes,6> CoccolithSimulation::supportedStokes({Stokes(1,1,0,0),Stokes(1,-1,0,0),Stokes(1,0,1,0),\
 Stokes(1,0,-1,0), Stokes(1,0,0,1),Stokes(1,0,0,-1)});
 Stokes CoccolithSimulation::incStoke;

double CoccolithSimulation::amplitudeScalingFactor = 1.0;

// Only for debugging
double sigmaTest( const meep::vec &r )
{
  return 1.0;
}

CoccolithSimulation::CoccolithSimulation(){};
CoccolithSimulation::~CoccolithSimulation()
{
  delete srcVol; srcVol=NULL;
  // Sources and flux are deleted in the destructor of field

  // Apperently it looks like field needs to be deleted before struc. Bug??
  delete field; field=NULL;
  delete struc; struc=NULL;
  delete dftVolTransmit; dftVolTransmit=NULL;
  delete dftVolRefl; dftVolRefl = NULL;
  delete dftVolBox; dftVolBox = NULL;
  delete transmitFlux; transmitFlux = NULL;
  delete reflFlux; reflFlux = NULL;
  delete fluxBox; fluxBox = NULL;
  delete file; file=NULL;
  delete sourceTime;sourceTime=NULL;

  delete monitor1; monitor1=NULL;
  delete monitor2; monitor2=NULL;
  delete faces; faces=NULL;
  delete n2fBox; n2fBox=NULL;
  #ifdef HAVE_LIB_VISA
    if ( plots != NULL ) delete plots; plots=NULL;
  #endif
}

cdouble CoccolithSimulation::amplitude( const meep::vec &r )
{
  cdouble im(0.0,1.0);
  if ( incStoke.U == -1 ) return -1.0*amplitudeScalingFactor;
  return 1.0*amplitudeScalingFactor; // In MEEP: & is dot product
}

cdouble CoccolithSimulation::amplitude2( const meep::vec &r )
{
  const double PI = acos(-1.0);
  cdouble im(0.0,1.0);
  if ( incStoke.V == 1 )
  {
    //return exp(im*PI/2.0)*amplitudeScalingFactor;
    return exp(-im*PI/2.0)*amplitudeScalingFactor;
  }
  else if ( incStoke.V == -1 )
  {
    //return exp(-im*PI/2.0)*amplitudeScalingFactor;
    return exp(+im*PI/2.0)*amplitudeScalingFactor;
  }
  return -1.0*amplitudeScalingFactor;
}

void CoccolithSimulation::setMainPropagationDirection( MainPropDirection_t propDir )
{
  propagationDir = propDir;
  addSourceVolume();
}

void CoccolithSimulation::addSourceVolume()
{
  switch( propagationDir )
  {
    case MainPropDirection_t::X:
      crn1 = meep::vec( getSrcPos(), gdvol.ymin(), gdvol.zmin() );
      crn2 = meep::vec( getSrcPos(), gdvol.ymax(), gdvol.zmax() );
      break;
    case MainPropDirection_t::Y:
      crn1 = meep::vec( gdvol.xmin(), getSrcPos(), gdvol.zmin() );
      crn2 = meep::vec( gdvol.xmax(), getSrcPos(), gdvol.zmax() );
      break;
    case MainPropDirection_t::Z:
      crn1 = meep::vec( gdvol.xmin(), gdvol.ymin(), getSrcPos() );
      crn2 = meep::vec( gdvol.xmax(), gdvol.ymax(), getSrcPos() );
      break;
  }

  if ( srcVol != NULL ) delete srcVol;
  srcVol = new meep::volume( crn1, crn2 );
}

void CoccolithSimulation::addStructure()
{

  if ( sourceTime == NULL )
  {
    throw( runtime_error("The source profile needs to be initialized before the structure is added!\n") );
  }
  if ( struc != NULL )
  {
    delete struc;
  }


  struc = new meep::structure( gdvol, *material, meep::pml( getPMLThickness() ) );
  struc->set_output_directory( outdir.c_str() );
  //struc->add_susceptibility( sigmaTest, meep::E_stuff, meep::lorentzian_susceptibility(0.1, 0.0) );

  if ( sellmeier != NULL )
  {
    updateStructure();
  }
}

double CoccolithSimulation::getWavelength() const
{
  return 1.0/centerFrequency; // In MEEP units c=1
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
  if ( geoIsInitialized )
  {
    throw( runtime_error("Geometry is already initialized! PML layers may behave strangely if the frequencies change!") );
  }
  centerFrequency = freq;
  freqWidth = fwidth;

  if ( source != NULL ) delete source;
  source = new meep::gaussian_src_time( freq, fwidth );
  sourceTime = source; // Default is to use a gaussian source
}

void CoccolithSimulation::addSource()
{
  if ( field == NULL )
  {
    throw( runtime_error("Field needs to be added before a source is added!\n") );
  }
  if ( sourceTime == NULL )
  {
    throw( runtime_error("A source profile needs to be added before a source is added! Call initSource!\n") );
  }
  if ( srcVol == NULL )
  {
    throw( runtime_error("A source volume needs to be added before a source is added!\n") );
  }

  switch( propagationDir )
  {
    case MainPropDirection_t::X:
      fieldComp = meep::component::Ey;
      secondComp = meep::Ez;
      break;
    case MainPropDirection_t::Y:
      fieldComp = meep::component::Ez;
      secondComp = meep::component::Ex;
      break;
    case MainPropDirection_t::Z:
      fieldComp = meep::component::Ex;
      secondComp = meep::component::Ey;
      break;
  }

  bool requireTwoSources = ( incStoke.U != 0 ) || ( incStoke.V != 0 );
  if ( initialStokesVectorSet && requireTwoSources )
  {
    meep::master_printf("Using two sources to obtain polarization\n");
    field->add_volume_source( secondComp, *sourceTime, *srcVol, amplitude2 );
    amplitudeScalingFactor = 1.0/sqrt(2.0);
  }

  if (( incStoke.Q == 1 ) || ( incStoke.Q == 0 ))
  {
    field->add_volume_source( fieldComp, *sourceTime, *srcVol, amplitude );
  }
  else
  {
    field->add_volume_source( secondComp, *sourceTime, *srcVol, amplitude );
  }
}

void CoccolithSimulation::init()
{
  if ( material == NULL )
  {
    throw( runtime_error("No material set! Needs to be set before call init()!") );
  }
  domainInfo();
  initializeGeometry();

  material->setDomainSize( gdvol, getPMLThickness()+additionalVaccumLayerPx );
  addSourceVolume();
  addStructure();
  addFields();
  addSource();
  addFluxPlanes();
  setMonitorPlanes();
  isInitialized = true;
}

double CoccolithSimulation::getPMLThickness() const
{
  return pmlThicknessInWavelengths*getWavelength();
}

void CoccolithSimulation::addFluxPlanes()
{

  if ( sourceTime == NULL )
  {
    throw( runtime_error("A source profile needs to be set before flux planes are added!") );
  }

  meep::vec transDFTCrn1;
  meep::vec transDFTCrn2;
  meep::vec reflDFTCrn1;
  meep::vec reflDFTCrn2;

  // Coordinates for flux box
  unsigned int lowerCrn[3];
  unsigned int upperCrn[3];

  // Shift the computed quantities to take into account the PML layers
  for ( unsigned int i=0;i<3;i++ )
  {
    lowerCrn[i] += getPMLThickness();
    upperCrn[i] += getPMLThickness();
  }
  //material->boundingBox( lowerCrn, upperCrn );
  //meep::vec boxCrn1( lowerCrn[0]-1.0, lowerCrn[1]-1.0, lowerCrn[2]-1.0 );
  //meep::vec boxCrn2( upperCrn[0]+1.0, upperCrn[1]+1.0, upperCrn[2]+1.0 );

  meep::vec boxCrn1;
  meep::vec boxCrn2;

  switch ( propagationDir )
  {
    case MainPropDirection_t::X:
      transDFTCrn1 = meep::vec( getTransFluxPos(), gdvol.ymin(), gdvol.zmin() );
      transDFTCrn2 = meep::vec( getTransFluxPos(), gdvol.ymax(), gdvol.zmax() );
      reflDFTCrn1 = meep::vec( getSrcFluxPos(), gdvol.ymin(), gdvol.zmin() );
      reflDFTCrn2 = meep::vec( getSrcFluxPos(), gdvol.ymax(), gdvol.zmax() );
      boxCrn1 = meep::vec( getTransFluxPos(), gdvol.ymin()+getPMLThickness(), gdvol.zmin()+getPMLThickness() );
      boxCrn2 = meep::vec( getSrcFluxPos(), gdvol.ymax()-getPMLThickness(), gdvol.zmax()-getPMLThickness() );
      break;
    case MainPropDirection_t::Y:
      transDFTCrn1 = meep::vec( gdvol.xmin(), getTransFluxPos(), gdvol.zmin() );
      transDFTCrn2 = meep::vec( gdvol.xmax(), getTransFluxPos(), gdvol.zmax() );
      reflDFTCrn1 = meep::vec( gdvol.xmin(), getSrcFluxPos(), gdvol.zmin() );
      reflDFTCrn2 = meep::vec( gdvol.xmax(), getSrcFluxPos(), gdvol.zmax() );
      boxCrn1 = meep::vec( gdvol.xmin()+getPMLThickness(), getTransFluxPos(), gdvol.zmin()+getPMLThickness() );
      boxCrn2 = meep::vec( gdvol.xmax()-getPMLThickness(), getTransFluxPos(), gdvol.zmax()-getPMLThickness() );
      break;
    case MainPropDirection_t::Z:
      transDFTCrn1 = meep::vec( gdvol.xmin(), gdvol.ymin(), getTransFluxPos() );
      transDFTCrn2 = meep::vec( gdvol.xmax(), gdvol.ymax(), getTransFluxPos() );
      reflDFTCrn1 = meep::vec( gdvol.xmin(), gdvol.ymin(), getSrcFluxPos() );
      reflDFTCrn2 = meep::vec( gdvol.xmax(), gdvol.ymax(), getSrcFluxPos() );
      boxCrn1 = meep::vec( gdvol.xmin()+getPMLThickness(), gdvol.ymin()+getPMLThickness(),  getSrcFluxPos() );
      boxCrn2 = meep::vec( gdvol.xmax()-getPMLThickness(), gdvol.ymax()-getPMLThickness(), getTransFluxPos() );
      break;
  }

  if ( dftVolTransmit != NULL ) delete dftVolTransmit;
  dftVolTransmit = new meep::volume( transDFTCrn1, transDFTCrn2 );

  transmitFlux = new meep::dft_flux( field->add_dft_flux_plane( *dftVolTransmit, centerFrequency-freqWidth/2.0, centerFrequency+freqWidth/2.0, nfreq ) );

  dftVolRefl = new meep::volume( reflDFTCrn1, reflDFTCrn2 );
  dftVolBox = new meep::volume( boxCrn1, boxCrn2 );
  reflFlux = new meep::dft_flux( field->add_dft_flux_plane( *dftVolRefl, centerFrequency-freqWidth/2.0, centerFrequency+freqWidth/2.0, nfreq ) );
  fluxBox = new meep::dft_flux( field->add_dft_flux_box( *dftVolBox, centerFrequency-freqWidth/2.0, centerFrequency+freqWidth/2.0, nfreq ) );

  if ( computeAsymmetryFactor ) addN2FPlanes( *dftVolBox );

  if ( !material->isReferenceRun() )
  {
    // There should be a backup file with the fluxes stored from the reference run
    reflFlux->load_hdf5( *field, reflFluxPlaneBackup.c_str() );
    fluxBox->load_hdf5( *field, reflFluxBoxBackup.c_str() );
    if ( n2fBox != NULL ) n2fBox->load_hdf5( *field, n2fBoxBackup.c_str() );
    reflFlux->scale_dfts(-1.0);
    fluxBox->scale_dfts(-1.0);
    if ( n2fBox != NULL ) n2fBox->scale_dfts(-1.0);
  }
}


double CoccolithSimulation::getLowerBorderInPropDir() const
{
  switch( propagationDir )
  {
    case MainPropDirection_t::X:
      return gdvol.xmin() + pmlThicknessInWavelengths*getWavelength();
    case MainPropDirection_t::Y:
      return gdvol.ymin() + pmlThicknessInWavelengths*getWavelength();
    case MainPropDirection_t::Z:
      return gdvol.zmin() + pmlThicknessInWavelengths*getWavelength();
  }
}

double CoccolithSimulation::getUpperBorderInPropDir() const
{
  switch( propagationDir )
  {
    case MainPropDirection_t::X:
      return gdvol.xmax() - pmlThicknessInWavelengths*getWavelength();
    case MainPropDirection_t::Y:
      return gdvol.ymax() - pmlThicknessInWavelengths*getWavelength();
    case MainPropDirection_t::Z:
      return gdvol.zmax() - pmlThicknessInWavelengths*getWavelength();
  }
}
double CoccolithSimulation::getSrcPos() const
{
  switch( srcPos )
  {
    case SourcePosition_t::TOP:
      return getLowerBorderInPropDir();
    case SourcePosition_t::BOTTOM:
      // In this case light propagates in the negative direction
      return getUpperBorderInPropDir();
  }
}

double CoccolithSimulation::getSrcFluxPos() const
{
  if ( material == NULL )
  {
    throw( runtime_error("Material needs to be set before call getSrcFluxPos()!") );
  }
  switch( srcPos )
  {
    case SourcePosition_t::TOP:
      return getSrcPos() + 5.0;//*material->getVoxelSize();
    case SourcePosition_t::BOTTOM:
      return getSrcPos() - 5.0;//*material->getVoxelSize();
  }
}

void CoccolithSimulation::setMonitorPlanes()
{
  if ( material == NULL )
  {
    throw( runtime_error("Material needs to be set before calling setMonitorPlanes()!") );
  }
  if ( monitor1 != NULL ) delete monitor1;
  if ( monitor2 != NULL ) delete monitor2;

  double dx = static_cast<double>(material->sizeX())/static_cast<double>(nMonitorX);
  double dy = static_cast<double>(material->sizeY())/static_cast<double>(nMonitorY);
  double dz = static_cast<double>(material->sizeZ())/static_cast<double>(nMonitorZ);

  meep::vec vx = meep::vec(dx, 0.0, 0.0 );
  meep::vec vy = meep::vec(0.0, dy, 0.0 );
  meep::vec vz = meep::vec(0.0, 0.0, dz );

  if ( meep::am_master() )
  {
    clog << "Setting up monitors and computing projections of epsilon...\n";
  }

  switch ( propagationDir )
  {
    case MainPropDirection_t::X:
      monitor1 = new FieldMonitor( nMonitorY, nMonitorX );
      monitor2 = new FieldMonitor( nMonitorZ, nMonitorX );
      monitor1->setName( "XY-plane" );
      monitor2->setName( "XZ-plane" );
      monitor1->setDisplacementVectors( vy, vx );
      monitor2->setDisplacementVectors( vz, vx );
      bkg1.set_size( nMonitorY, nMonitorX );
      bkg2.set_size( nMonitorZ, nMonitorX );
      projectedEpsilon( bkg1, IntegrationDir_t::Z );
      projectedEpsilon( bkg2, IntegrationDir_t::Y );
      break;
    case MainPropDirection_t::Y:
      monitor1 = new FieldMonitor( nMonitorY, nMonitorX );
      monitor2 = new FieldMonitor( nMonitorZ, nMonitorY );
      monitor1->setName( "XY-plane" );
      monitor2->setName( "YZ-plane" );
      monitor1->setDisplacementVectors( vy, vx );
      monitor2->setDisplacementVectors( vz, vy );
      bkg1.set_size( nMonitorY, nMonitorX );
      bkg2.set_size( nMonitorZ, nMonitorY );
      projectedEpsilon( bkg1, IntegrationDir_t::Z );
      projectedEpsilon( bkg2, IntegrationDir_t::X );
      break;
    case MainPropDirection_t::Z:
      monitor1 = new FieldMonitor( nMonitorZ, nMonitorX );
      monitor2 = new FieldMonitor( nMonitorZ, nMonitorY );
      monitor1->setName( "XZ-plane" );
      monitor2->setName( "YZ-plane" );
      material->projectionXZ( bkg1 );
      material->projectionYZ( bkg2 );
      monitor1->setDisplacementVectors( vz, vx );
      monitor2->setDisplacementVectors( vz, vy );
      bkg1.set_size( nMonitorZ, nMonitorX );
      bkg2.set_size( nMonitorZ, nMonitorY );
      projectedEpsilon( bkg1, IntegrationDir_t::Y );
      projectedEpsilon( bkg2, IntegrationDir_t::X );
      break;
  }

  monitor1->setOrigin( gdvol.center() );
  monitor2->setOrigin( gdvol.center() );

  if ( realTimeVisualization )
  {
    #ifdef HAVE_LIB_VISA
      if ( plots != NULL ) delete plots;

      plots = new visa::WindowHandler();
      // Add plots
      plots->addPlot( monitor1->getName().c_str() );
      plots->addPlot( monitor2->getName().c_str() );
      plots->setLayout( 1, 2 );
      plots->useSeparateDrawing(); // Do not show result on screen before show
    #endif
  }

  if ( meep::am_master() )
  {
    clog << "Finished\n";
  }

}

void CoccolithSimulation::run()
{
  if ( !isInitialized )
  {
    throw( runtime_error("Simulation is not initialized!") );
  }

  unsigned int iter = 0;
  double simStop = source->last_time() + estimatedTimeToPropagateAcrossDomain();

  if ( userOverridedEndTime )
  {
    simStop = tEnd;
  }

  if ( meep::am_master() )
  {
    clog << "End time: " << simStop << endl;
  }

  while ( field->time() < simStop )
  {
    field->step();
    if (( iter%plotUpdateFreq == 0 ) && ( realTimeVisualization ))
    {
      visualize();
    }
    iter++;
  }
}

void CoccolithSimulation::visualize()
{
  #ifdef HAVE_LIB_VISA
    if (( monitor1 == NULL ) || ( monitor2 == NULL ))
    {
      throw (runtime_error("No monitors set for visulization!") );
    }

    assert( plots != NULL );

    monitor1->setIntensity( *field );
    monitor2->setIntensity( *field );
    // Set colorbar limits

    typedef visa::Colormaps::Colormap_t cmap_t;

    visa::Visualizer& plt1 = plots->get( monitor1->getName().c_str() );
    visa::Visualizer& plt2 = plots->get( monitor2->getName().c_str() );

    plt1.setColorLim( monitor1->get().min(), monitor1->get().max() );
    plt2.setColorLim( monitor2->get().min(), monitor2->get().max() );
    plt1.setOpacity( 1.0 );
    plt2.setOpacity( 1.0 );
    plt1.setCmap( cmap_t::NIPY_SPECTRAL );
    plt2.setCmap( cmap_t::NIPY_SPECTRAL );
    plt1.setImg( monitor1->get() );
    plt2.setImg( monitor2->get() );
    plots->draw();

    // Overlay refractive index profile
    plt1.setCmap( cmap_t::GREYSCALE );
    plt2.setCmap( cmap_t::GREYSCALE );
    plt1.setColorLim( 0.99*bkg1.min(), 1.01*bkg1.max() );
    plt2.setColorLim( 0.99*bkg2.min(), 1.01*bkg2.max() );
    plt1.setOpacity( 0.6 );
    plt2.setOpacity( 0.6 );
    plt1.setImg( bkg1 );
    plt2.setImg( bkg2 );
    plots->draw();

    plots->show();
  #endif
}

void CoccolithSimulation::domainInfo() const
{
  if ( meep::am_master() )
  {
    if ( material == NULL )
    {
      throw( runtime_error("Material needs to set be before calling domainInfo()!") );
    }
    double speedOfLight = 2.997E8; // nm/ns
    double L = material->getVoxelSize(); // In nm
    double PI = acos(-1.0);
    cout << "-----------------------------------------------------------\n";
    cout << "Domain info:\n";
    cout << "SizeX: " << material->sizeX()*L << " nm, SizeY: " << material->sizeY()*L << " nm. SizeZ: " << material->sizeZ()*L << " nm\n";
    cout << "dx=dy=dz = " << L/resolution << " nm\n";
    cout << "PML thickness: " << getPMLThickness()*L << " nm\n";
    cout << "Main frequency: " << centerFrequency*speedOfLight/(L*1000.0) << " (THz)\n";
    cout << "Frequency width: " << freqWidth*speedOfLight/(L*1000.0) << " (THz)\n";
    cout << "Wavelength: " << getWavelength()*L << " nm\n";
    cout << "Compute asymmetry factor: " << computeAsymmetryFactor << endl;
    cout << "Compute Stokes parameters: " << computeStokesParameters << endl;
  }
}

void CoccolithSimulation::projectedEpsilon( arma::mat &values, IntegrationDir_t dir )
{
  if ( material == NULL )
  {
    throw( runtime_error("Material needs to be set before calling projectedEpsilon()!") );
  }
  double dx, dy, dz;
  unsigned int nIntegr = 0;
  Plane_t plane;
  switch( dir )
  {
    case IntegrationDir_t::X:
      plane = Plane_t::YZ;
      nIntegr = nMonitorX;
      dx = ( gdvol.xmax()-gdvol.xmin() )/nMonitorX;
      dy = ( gdvol.ymax() - gdvol.ymin() )/values.n_cols;
      dz = ( gdvol.zmax() - gdvol.zmin() )/values.n_rows;
      break;
    case IntegrationDir_t::Y:
      nIntegr = nMonitorY;
      plane = Plane_t::XZ;
      dx = ( gdvol.xmax() - gdvol.xmin() )/values.n_cols;
      dy = ( gdvol.ymax() - gdvol.ymin() )/nMonitorY;
      dz = ( gdvol.zmax() - gdvol.zmin() )/values.n_rows;
      break;
    case IntegrationDir_t::Z:
      plane = Plane_t::XY;
      nIntegr = nMonitorZ;
      dx = ( gdvol.xmax() - gdvol.xmin() )/values.n_cols;
      dy = ( gdvol.ymax() - gdvol.ymin() )/values.n_rows;
      dz = ( gdvol.zmax() - gdvol.zmin() )/nMonitorZ;
      break;
  }

  values.fill(0.0);

  #pragma omp parallel
  {
    meep::vec pos;
    #pragma omp for
    for ( unsigned int i=0;i<values.n_cols*values.n_rows;i++ )
    {
      unsigned int row = i%values.n_rows;
      unsigned int col = i/values.n_rows;
      for ( unsigned int k=0;k<nIntegr;k++ )
      {
        getPos( col, row, plane, dx, dy, dz, pos );
        switch( dir )
        {
          case IntegrationDir_t::X:
            pos += meep::vec( gdvol.xmin() + k*dx, 0.0, 0.0 );
            break;
          case IntegrationDir_t::Y:
            pos += meep::vec( 0.0, gdvol.ymin() + k*dy, 0.0 );
            break;
          case IntegrationDir_t::Z:
            pos += meep::vec( 0.0, 0.0, gdvol.zmin() + k*dz );
            break;
        }
        values(col,row) += material->eps( pos );
      }
      values(col,row) /= nIntegr;
    }
  }
}

void CoccolithSimulation::getPos( unsigned int row, unsigned int col, Plane_t plane, double dx, double dy, double dz, meep::vec &res ) const
{
  switch( plane )
  {
    case Plane_t::XY:
      // col: x, row: y
      res = meep::vec( gdvol.xmin() + col*dx, gdvol.ymin() + row*dy,0.0);
      return;
    case Plane_t::XZ:
      // col: x, row: y
      res = meep::vec( gdvol.xmin() + col*dx, 0.0, gdvol.zmin() + row*dz );
      return;
    case Plane_t::YZ:
      // col: y, row: z
      res = meep::vec( 0.0, gdvol.ymin() + col*dy, gdvol.zmin() + row*dz );
      return;
  }
}

void CoccolithSimulation::setUID()
{
  if ( uid == "" )
  {
    time_t t = time(NULL);
    auto tm = *localtime(&t);
    stringstream ss;
    // Skip seconds as different processes reach this point slightly different, but a minute change in minutes is more unlikely
    ss << put_time( &tm, "%Y%m%d_%H%M");
    uid = ss.str();
  }
}
void CoccolithSimulation::exportResults()
{
  if ( material == NULL )
  {
    throw( runtime_error("Material needs to be set before calling exportResults()!") );
  }

  setUID();
  stringstream ss;
  if ( prefix == "" )
  {
    prefix = "defaultFilename";
  }
  ss << prefix << "_" << uid; // File extension added automatically

  // Send the filename to the other processes
  int root = 0;
  int size = ss.str().size();
  char buffer[size];
  memcpy(buffer, ss.str().c_str(),size);
  meep::broadcast(root, buffer, size );
  if ( file == NULL )
  {
    if ( material->isReferenceRun() )
    {
      file = field->open_h5file( buffer );
      saveGeometry();
      saveDFTParameters();
    }
    else
    {
      file = field->open_h5file( buffer, meep::h5file::READWRITE );
    }
  }
  if ( !material->isReferenceRun() )
  {
    field->output_hdf5( meep::Dielectric, gdvol.surroundings(), file, false, true );
    file->prevent_deadlock();

    field->output_hdf5( fieldComp, gdvol.surroundings(), file, false, true );
    file->prevent_deadlock();

    field->output_hdf5( secondComp, gdvol.surroundings(), file, false, true );
    file->prevent_deadlock();
    farFieldQuantities();
  }
  saveDFTSpectrum();
  //saveDFTStokes();

  if (( meep::am_master() ) && (!material->isReferenceRun() ))
  {
    clog << "Results written to " << ss.str() << endl;
  }
}

void CoccolithSimulation::farFieldQuantities()
{
  if ( file == NULL )
  {
    file = field->open_h5file( "data/farfield" );
  }
  if ( computeAsymmetryFactor )
  {
    double PI = acos(-1.0);
    double nFreq = n2fBox->Nfreq;
    double ffFreq[3] = {n2fBox->freq_min/(2.0*PI), n2fBox->dfreq/(2.0*PI), nFreq};

    int nFFfreq = 3;
    file->write("ffFreq", 1, &nFFfreq, ffFreq, false );
    file->prevent_deadlock();

    vector<double> asym;
    scatteringAssymmetryFactor( asym, 40000.0, gaussLegendreOrder );
    int length = asym.size();

    file->write( "Asymmetry", 1, &length, &asym[0], false );
    file->prevent_deadlock();

    int gLgorder = gaussLegendreOrder;
    int nfreqCast = radialPoyntingVector.n_cols;
    int ntheta = radialPoyntingVector.n_rows;
    int rank[2] = {ntheta, nfreqCast};
    arma::mat radPoyntingTrans = radialPoyntingVector.t(); // Is MEEP changing the memory layout?
    file->write( "Sr", 2, rank, radPoyntingTrans.memptr(), false );
    file->prevent_deadlock();

    file->write( "SrTheta", 1, &gLgorder, &thetaValues[0], false );
    file->prevent_deadlock();

    int size = phiValues.size();
    file->write("phi", 1, &size, &phiValues[0], false );
    file->prevent_deadlock();

    if ( computeStokesParameters )
    {
      rank[0] = stokesI[0].n_rows;
      rank[1] = stokesI[0].n_cols;
      bool saveStokesParameters = (stokesI[0].n_rows == stokesQ[0].n_rows ) && ( stokesI[0].n_rows == stokesU[0].n_rows ) && \
                                  (stokesI[0].n_rows == stokesV[0].n_rows ) && ( stokesI[0].n_cols == stokesQ[0].n_cols ) && \
                                  (stokesI[0].n_cols == stokesU[0].n_cols ) && ( stokesI[0].n_cols == stokesV[0].n_cols );

      if ( saveStokesParameters )
      {
        rank[0] = stokesIAzim.n_rows;
        rank[1] = stokesIAzim.n_cols;
        arma::mat stokes = stokesIAzim.t();
        file->write("stokesIAvg", 2, rank, stokes.memptr(), false );
        file->prevent_deadlock();

        stokes = stokesQAzim.t();
        file->write("stokesQAvg", 2, rank, stokes.memptr(), false );
        file->prevent_deadlock();

        stokes = stokesUAzim.t();
        file->write("stokesUAvg", 2, rank, stokes.memptr(), false );
        file->prevent_deadlock();

        stokes = stokesVAzim.t();
        file->write("stokesVAvg", 2, rank, stokes.memptr(), false );
        file->prevent_deadlock();
        //n2fBox->save_hdf5( *field, "data/nearToFarScattering" );

        saveStokesPhiTheta();

      }
      else
      {
        meep::master_printf("Warning! The dimensions of the different Stokes parameters does not match!\n");
      }
    }
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

void CoccolithSimulation::saveDFTSpectrum()
{
  //assert( file != NULL );
  assert( material != NULL );

  string newpref("");
  if ( material->isReferenceRun() )
  {
    int length = transmitFlux->Nfreq;
    file->write( "spectrumReference", 1, &length, transmitFlux->flux(), false );
    file->prevent_deadlock();
    file->write( "boxFluxRef", 1, &length, fluxBox->flux(), false );
    file->prevent_deadlock();
    file->write("reflPlaneRef", 1, &length, reflFlux->flux(), false );
    file->prevent_deadlock();
    // Save backup fields
    fluxBox->save_hdf5( *field, reflFluxBoxBackup.c_str() );
    reflFlux->save_hdf5( *field, reflFluxPlaneBackup.c_str() );
    n2fBox->save_hdf5( *field, n2fBoxBackup.c_str() );
  }
  else
  {
    int length = transmitFlux->Nfreq;
    file->write( "spectrumTransmitted", 1, &length, transmitFlux->flux(), false );
    file->prevent_deadlock();
    file->write( "boxFluxScat", 1, &length, fluxBox->flux(), false );
    file->prevent_deadlock();
    file->write("reflPlaneScat", 1, &length, reflFlux->flux(), false );
    file->prevent_deadlock();
  }
}

void CoccolithSimulation::saveDFTParameters()
{
  const int nParams = 3;
  double stat[nParams];
  stat[0] = transmitFlux->freq_min;
  stat[1] = transmitFlux->dfreq;
  stat[2] = transmitFlux->Nfreq;

  // Write parameters
  file->write( "spectrumFreqs", 1, &nParams, stat, false );
  file->prevent_deadlock();

  // Write geometrical positions
  meep::vec mincrn = dftVolTransmit->get_min_corner();
  meep::vec maxcrn = dftVolTransmit->get_max_corner();
  int nCrns = 6;
  double cornersTrans[nCrns] = {mincrn.x(), mincrn.y(), mincrn.z(), maxcrn.x(), maxcrn.y(), maxcrn.z()};
  file->write( "spectrumCornersTransmit", 1, &nCrns, cornersTrans, false );
  file->prevent_deadlock();
}

void CoccolithSimulation::saveGeometry()
{
  assert( file != NULL );

  meep::vec mincrn = srcVol->get_min_corner();
  meep::vec maxcrn = srcVol->get_max_corner();
  int nCrn = 6;
  double corners[6] = {mincrn.x(), mincrn.y(), mincrn.z(), maxcrn.x(), maxcrn.y(), maxcrn.z()};
  file->write( "geometrySourceVolume", 1, &nCrn, corners, false );
  file->prevent_deadlock();


  int ndim = 3;
  double incWaveVec[ndim] = {waveVec.x(), waveVec.y(), waveVec.z()};
  file->write( "geometryIncWavevector", 1, &ndim, incWaveVec, false );
  file->prevent_deadlock();

  int single = 1;
  double wavelength = getWavelength();
  file->write( "geometryWavelength", 1, &single, &wavelength, false );
  file->prevent_deadlock();


  int nFreqParam = 2;
  double freqParam[nFreqParam] = {centerFrequency, freqWidth};
  file->write( "geometryAngFreq", 1, &nFreqParam, freqParam, false );
  file->prevent_deadlock();


  double domainSize[nCrn] = {gdvol.xmin(), gdvol.xmax(), gdvol.ymin(), gdvol.ymax(), gdvol.zmin(), gdvol.zmax()};
  file->write( "geometryDomainSize", 1, &nCrn, domainSize, false );
  file->prevent_deadlock();


  double pmlThickness = getPMLThickness();
  file->write( "geometryPmlT", 1, &single, &pmlThickness, false );
  file->prevent_deadlock();

  double vxsize = material->getVoxelSize();
  file->write( "vxSizeNM", 1, &single, &vxsize, false );
  file->prevent_deadlock();

  int stkSize = 4;
  double incidentStokesVector[stkSize];
  incidentStokesVector[0] = incStoke.I;
  incidentStokesVector[1] = incStoke.Q;
  incidentStokesVector[2] = incStoke.U;
  incidentStokesVector[3] = incStoke.V;
  file->write("incStokes", 1, &stkSize, incidentStokesVector, false);
  file->prevent_deadlock();

  if ( dftVolBox != NULL )
  {
    meep::vec bCrn1 = dftVolBox->get_min_corner();
    meep::vec bCrn2 = dftVolBox->get_max_corner();
    double bxCrn[nCrn] = {bCrn1.x(),bCrn1.y(),bCrn1.z(),bCrn2.x(),bCrn2.y(),bCrn2.z()};
    file->write("boxGeo", 1,&nCrn,bxCrn,false);
    file->prevent_deadlock();
  }

  if ( dftVolRefl != NULL )
  {
    meep::vec refPCrn1 = dftVolRefl->get_min_corner();
    meep::vec refPCrn2 = dftVolRefl->get_max_corner();
    double refVolCrn[nCrn] = {refPCrn1.x(),refPCrn1.y(),refPCrn1.z(),refPCrn2.x(),refPCrn2.y(),refPCrn2.z()};
    file->write("refPGeo", 1, &nCrn, refVolCrn, false );
    file->prevent_deadlock();
  }

  if ( dftVolTransmit != NULL )
  {
    meep::vec trPCrn1 = dftVolTransmit->get_min_corner();
    meep::vec trPCrn2 = dftVolTransmit->get_max_corner();
    double trVolCrn[nCrn] = {trPCrn1.x(),trPCrn1.y(),trPCrn1.z(),trPCrn2.x(),trPCrn2.y(),trPCrn2.z()};
    file->write("trPGeo", 1, &nCrn, trVolCrn, false );
    file->prevent_deadlock();
  }
}

void CoccolithSimulation::runWithoutScatterer()
{
  if ( material == NULL )
  {
    throw( runtime_error("Material needs to be set before calling runWithoutScatterer()!") );
  }
  material->setReferenceRun( true );
}

void CoccolithSimulation::runWithScatterer()
{
  if ( material == NULL )
  {
    throw( runtime_error("Material needs to be set before calling runWithScatterer()!") );
  }
  material->setReferenceRun( false );
}

void CoccolithSimulation::reset()
{
  field->reset();
  init();
}

void CoccolithSimulation::setEndTime( double newtime )
{
  tEnd = newtime;
  userOverridedEndTime = true;
}

double CoccolithSimulation::estimatedTimeToPropagateAcrossDomain() const
{
  if ( material == NULL )
  {
    throw( runtime_error("Material needs to be set before calling estimatedTimeToPropagateAcrossDomain()!") );
  }
  arma::uvec sizes(3);
  sizes(0) = material->sizeX();
  sizes(1) = material->sizeY();
  sizes(2) = material->sizeZ();
  return sizes.max();
}

void CoccolithSimulation::initializeGeometry()
{
  if ( material == NULL )
  {
    throw ( runtime_error("Material needs to be loaded before the geometry is initialize!") );
  }

  double pmlThickPx = getPMLThickness();
  double extra = 2.0*pmlThickPx + 2.0*additionalVaccumLayerPx;
  gdvol = meep::vol3d( material->sizeX()+extra, material->sizeY()+extra, material->sizeZ()+extra, resolution );
  geoIsInitialized = true;
}

void CoccolithSimulation::setPMLInWavelengths( double pmlInWav )
{
  if ( geoIsInitialized )
  {
    throw( runtime_error("Geometry is already initialized. It has no effect to change the PML thickness now!"));
  }
  pmlThicknessInWavelengths = pmlInWav;
}

void CoccolithSimulation::updateStructure()
{
  if ( struc == NULL )
  {
    throw (runtime_error("Structure is not initialized!") );
  }

  if ( sellmeier == NULL )
  {
    throw( runtime_error("Sellmeier material not specified!") );
  }
  const double PI = acos(-1.0);
  for ( unsigned int i=0;i<sellmeier->nLorentzians();i++ )
  {
    double sigma, omega0;
    sellmeier->getMEEPLorentzian( material->getVoxelSize()*1E-3, i, sigma, omega0 );
    VoxelSusceptibility matFunc( sigma, 0.0 );
    struc->add_susceptibility( matFunc, meep::E_stuff, meep::lorentzian_susceptibility(omega0/(2.0*PI), 0.0) );
  }

  if ( meep::am_master() )
  {
    clog << "Added " << sellmeier->nLorentzians() << " lorentzian susceptibilities to the MEEP structure...\n";
  }
}

void CoccolithSimulation::addN2FPlanes( const meep::volume &box )
{
  delete faces;
  delete n2fBox;
  faces = NULL;
  n2fBox = NULL;
  const meep::vec mincrn = box.get_min_corner();
  const meep::vec maxcrn = box.get_max_corner();

  // Face to store flux in negative x-direction
  meep::vec crn1 = mincrn;
  meep::vec crn2( mincrn.x(), maxcrn.y(), maxcrn.z() );
  faces = new meep::volume_list( meep::volume(crn1, crn2), meep::Sx, -1, faces );
  #ifdef DEBUG_N2F_INITIALIZATION
    meep::master_printf( "mincrn:%.1f,%.1f,%.1f, maxcrn=%.1f,%.1f,%.1f\n", crn1.x(), crn1.y(), crn1.z(), crn2.x(),crn2.y(),crn2.z() );
  #endif

  // Face to store flux in positive x-direction
  crn1 = meep::vec( maxcrn.x(), mincrn.y(), mincrn.z() );
  crn2 = maxcrn;
  faces = new meep::volume_list( meep::volume(crn1, crn2), meep::Sx, 1 , faces);
  #ifdef DEBUG_N2F_INITIALIZATION
    meep::master_printf( "mincrn:%.1f,%.1f,%.1f, maxcrn=%.1f,%.1f,%.1f\n", crn1.x(), crn1.y(), crn1.z(), crn2.x(),crn2.y(),crn2.z() );
  #endif

  // Face to store flux in negative y-direction
  crn1 = mincrn;
  crn2 = meep::vec( maxcrn.x(), mincrn.y(), maxcrn.z() );
  faces = new meep::volume_list( meep::volume(crn1, crn2), meep::Sy, -1, faces );
  #ifdef DEBUG_N2F_INITIALIZATION
    meep::master_printf( "mincrn:%.1f,%.1f,%.1f, maxcrn=%.1f,%.1f,%.1f\n", crn1.x(), crn1.y(), crn1.z(), crn2.x(),crn2.y(),crn2.z() );
  #endif

  // Face to store in positive y-direction
  crn1 = meep::vec( mincrn.x(), maxcrn.y(), mincrn.z() );
  crn2 = maxcrn;
  faces = new meep::volume_list( meep::volume(crn1, crn2), meep::Sy, 1 , faces);
  #ifdef DEBUG_N2F_INITIALIZATION
    meep::master_printf( "mincrn:%.1f,%.1f,%.1f, maxcrn=%.1f,%.1f,%.1f\n", crn1.x(), crn1.y(), crn1.z(), crn2.x(),crn2.y(),crn2.z() );
  #endif

  // Face to store flux in negative z-direction
  crn1 = mincrn;
  crn2 = meep::vec( maxcrn.x(), maxcrn.y(), mincrn.z() );
  faces = new meep::volume_list( meep::volume(crn1, crn2), meep::Sz, -1, faces );
  #ifdef DEBUG_N2F_INITIALIZATION
    meep::master_printf( "mincrn:%.1f,%.1f,%.1f, maxcrn=%.1f,%.1f,%.1f\n", crn1.x(), crn1.y(), crn1.z(), crn2.x(),crn2.y(),crn2.z() );
  #endif

  // Face to store flux in positive z-direction
  crn1 = meep::vec( mincrn.x(), mincrn.y(), maxcrn.z() );
  crn2 = maxcrn;
  faces = new meep::volume_list( meep::volume(crn1, crn2), meep::Sz, 1 , faces);
  #ifdef DEBUG_N2F_INITIALIZATION
    meep::master_printf( "mincrn:%.1f,%.1f,%.1f, maxcrn=%.1f,%.1f,%.1f\n", crn1.x(), crn1.y(), crn1.z(), crn2.x(),crn2.y(),crn2.z() );
  #endif

  n2fBox = new meep::dft_near2far( field->add_dft_near2far(faces, centerFrequency-freqWidth/2.0, centerFrequency+freqWidth/2.0, 5 ) );
}


void CoccolithSimulation::scatteringAssymmetryFactor( vector<double> &g, double R, unsigned int Nsteps )
{
  meep::master_printf( "Computing assymmetry factor...\n" );
  radialPoyntingVector.set_size( Nsteps, n2fBox->Nfreq );
  thetaValues.resize(Nsteps);
  g.resize( n2fBox->Nfreq );
  Ephi.set_size( Nsteps, numberOfAzimuthalSteps );
  Etheta.set_size( Nsteps, numberOfAzimuthalSteps );

  if ( computeStokesParameters )
  {
    stokesI.resize( numberOfAzimuthalSteps );
    stokesQ.resize( numberOfAzimuthalSteps );
    stokesU.resize( numberOfAzimuthalSteps );
    stokesV.resize( numberOfAzimuthalSteps );
    stokesIInc.resize( numberOfAzimuthalSteps );
    stokesQInc.resize( numberOfAzimuthalSteps );
    stokesUInc.resize( numberOfAzimuthalSteps );
    stokesVInc.resize( numberOfAzimuthalSteps );
    for ( unsigned int i=0;i<numberOfAzimuthalSteps;i++ )
    {
      stokesI[i].set_size(Nsteps, n2fBox->Nfreq);
      stokesQ[i].set_size(Nsteps, n2fBox->Nfreq);
      stokesU[i].set_size(Nsteps, n2fBox->Nfreq);
      stokesV[i].set_size(Nsteps, n2fBox->Nfreq);
    }
    stokesIAzim.set_size(Nsteps,n2fBox->Nfreq);
    stokesQAzim.set_size(Nsteps,n2fBox->Nfreq);
    stokesUAzim.set_size(Nsteps,n2fBox->Nfreq);
    stokesVAzim.set_size(Nsteps,n2fBox->Nfreq);

    stokesIAzim.fill(0.0);
    stokesQAzim.fill(0.0);
    stokesUAzim.fill(0.0);
    stokesVAzim.fill(0.0);
  }
  vector<double> normalization(g.size());
  double PI = acos(-1.0);
  gsl_integration_glfixed_table *gslTab = gsl_integration_glfixed_table_alloc(Nsteps);

  vector<double> azimInt;
  fill( g.begin(),g.end(),0.0);
  fill( normalization.begin(), normalization.end(), 0.0 );
  for ( unsigned int i=0;i<Nsteps;i++ )
  {
    currentTheta = i;
    double theta;
    double weight;
    gsl_integration_glfixed_point( 0.0, PI, i, &theta, &weight, gslTab );

    int percentageDone = static_cast<int>(i*100.0/Nsteps);
    meep::master_printf("%d done...\r", percentageDone);
    //azimuthalIntagration( theta, R, Nsteps, azimInt);
    azimuthalIntagration( theta, R, numberOfAzimuthalSteps, azimInt);
    double cosTheta = 0.0;
    if ( redefineTheta() )
    {
      thetaValues[i] = PI-theta;
      cosTheta = cos(PI-theta);
    }
    else
    {
      thetaValues[i] = theta;
      cosTheta = cos(theta);
    }
    for ( unsigned int f=0;f<n2fBox->Nfreq;f++ )
    {
      g[f] += azimInt[f]*cosTheta*weight*sin(theta);
      normalization[f] += azimInt[f]*weight*sin(theta);
      radialPoyntingVector(i,f) = azimInt[f];
    }
  }
  if ( meep::am_master() )
  {
    clog << endl;
  }

  // Normalize
  for ( unsigned int i=0;i<g.size();i++ )
  {
    g[i] /= normalization[i];
  }
  gsl_integration_glfixed_table_free( gslTab );
  meep::master_printf( "Finished computing assymmetry factor\n" );
}

void CoccolithSimulation::azimuthalIntagration( double theta, double R, unsigned int Nsteps, vector<double> &res )
{
  // Determine the weights and evaluation points for phi [0,2*pi]
  gsl_integration_glfixed_table *gslTab = gsl_integration_glfixed_table_alloc(Nsteps);

  if ( res.size() != n2fBox->Nfreq ) res.resize(n2fBox->Nfreq);
  fill(res.begin(),res.end(),0.0);
  cdouble *results = NULL;
  double RsinTheta = R*sin(theta);
  double RcosTheta = R*cos(theta);
  double add = 0.0;
  double PI = acos(-1.0);

  // Apply Gauss-Legendre integation scheme
  for ( unsigned int n=0;n<Nsteps;n++ )
  {
    double phi, weight;
    gsl_integration_glfixed_point(0.0, 2.0*PI, n, &phi, &weight, gslTab );

    if ( currentTheta == 0 ) phiValues.push_back(phi);

    double x = RsinTheta*cos(phi);
    double y = RsinTheta*sin(phi);
    double z = RcosTheta;
    permumteToFitPropDir(x,y,z);

    results = n2fBox->farfield( meep::vec(x,y,z) );
    for ( unsigned int i=0;i<n2fBox->Nfreq;i++ )
    {
      //add = pow(abs(results[6*i]),2) + pow(abs(results[6*i+1]),2) + pow(abs(results[6*i+2]),2);
      res[i] += weight*phaseFunctionContribution(results+6*i);
    }

    if ( computeStokesParameters )
    {
      LocalStokes locStokes;
      Stokes incTransformed;
      getLocalStokes( theta, phi, results, locStokes, incTransformed );
      for ( unsigned int i=0;i<n2fBox->Nfreq;i++ )
      {
        stokesI[n](currentTheta,i) = locStokes.I[i];
        stokesQ[n](currentTheta,i) = locStokes.Q[i];
        stokesU[n](currentTheta,i) = locStokes.U[i];
        stokesV[n](currentTheta,i) = locStokes.V[i];
        stokesIAzim(currentTheta,i) += locStokes.I[i]*weight;
        stokesQAzim(currentTheta,i) += locStokes.Q[i]*weight;
        stokesUAzim(currentTheta,i) += locStokes.U[i]*weight;
        stokesVAzim(currentTheta,i) += locStokes.V[i]*weight;
      }
      if ( currentTheta == 0 )
      {
        stokesIInc[n] = incTransformed.I;
        stokesQInc[n] = incTransformed.Q;
        stokesUInc[n] = incTransformed.U;
        stokesVInc[n] = incTransformed.V;
      }

      Ephi(currentTheta,n) = locStokes.Ephi;
      Etheta(currentTheta,n) = locStokes.Etheta;

      //updateStokesParameters( results, currentTheta, weight );
    }
    delete results; results=NULL;
  }
  gsl_integration_glfixed_table_free( gslTab );
}

template <class T>
void CoccolithSimulation::permumteToFitPropDir( T &x, T &y, T &z ) const
{
  T copyX = x;
  T copyY = y;
  T copyZ = z;
  switch( propagationDir )
  {
    case MainPropDirection_t::X:
      x = copyZ;
      y = copyX;
      z = copyY;
      return;
    case MainPropDirection_t::Y:
      x = copyY;
      y = copyZ;
      z = copyX;
      return;
    case MainPropDirection_t::Z:
      return;
  }
}

// Explicit initialization of allowed template arguments
template void CoccolithSimulation::permumteToFitPropDir<double>( double &x, double &y, double &z ) const;
template void CoccolithSimulation::permumteToFitPropDir<cdouble>( cdouble &x, cdouble &y, cdouble &z ) const;

void CoccolithSimulation::updateStokesParameters( const cdouble EH[], unsigned int evalPointIndx, double weight )
{
  /*
  cdouble E1, E2;
  for( unsigned int i=0;i<;i++ )
  {
    E1 = EH[6*i+meep::component_index(fieldComp)];
    E2 = EH[6*i+meep::component_index(secondComp)];
    stokesI(evalPointIndx, i) += (pow( abs(E1),2 )+pow( abs(E2),2) )*weight;
    stokesU(evalPointIndx, i) += (pow(abs(E1),2)-pow(abs(E2),2))*weight;
    stokesQ(evalPointIndx, i) += 2.0*(E1*conj(E2)).real()*weight;
    stokesV(evalPointIndx, i) -= 2.0*(E1*conj(E2)).imag()*weight;
  }*/
}

double CoccolithSimulation::phaseFunctionContribution( const cdouble EH[3] ) const
{
  double contribution = 0.0;
  switch ( propagationDir )
  {
    case MainPropDirection_t::X:
      return pow(abs(EH[1]),2) + pow(abs(EH[2]),2);
    case MainPropDirection_t::Y:
      return pow(abs(EH[0]),2) + pow(abs(EH[2]),2);
    case MainPropDirection_t::Z:
      return pow(abs(EH[1]),2) + pow(abs(EH[0]),2);
  }
}

bool CoccolithSimulation::redefineTheta() const
{
  // If the source is at the bottom, the pulse propagates along the negative axis
  // In this case theta=0 should correspond to propagation along the negative axis
  // Otherwise, theta=0 should correspond to propagation along the positive axis
  return srcPos == SourcePosition_t::BOTTOM;
}

void CoccolithSimulation::computeEvectorOrthogonalToPropagation( double theta, double phi, meep::vec &E1hat, meep::vec &E2hat ) const
{
  double x = sin(theta)*cos(phi);
  double y = sin(theta)*sin(phi);
  double z = cos(theta);
  permumteToFitPropDir(x,y,z);
  meep::vec khat( x, y, z );

  // Unit vector for theta hat
  x = cos(theta)*cos(phi);
  y = cos(theta)*sin(phi);
  z = -sin(theta);
  permumteToFitPropDir(x,y,z);
  E1hat = meep::vec(x,y,z);

  // Unit vector in phi hat
  x = -sin(phi);
  y = cos(phi);
  z = 0.0;
  permumteToFitPropDir(x,y,z);
  E2hat = meep::vec(x,y,z);

  // Some debug assertions, verify that the vectors are orthonormal
  double E1len = sqrt( E1hat.x()*E1hat.x() + E1hat.y()*E1hat.y() + E1hat.z()*E1hat.z() );
  double E2len = sqrt( E2hat.x()*E2hat.x() + E2hat.y()*E2hat.y() + E2hat.z()*E2hat.z() );
  double khatlen = sqrt( khat.x()*khat.x() + khat.y()*khat.y() + khat.z()*khat.z() );
  //meep::master_printf("%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n", abs(khat&E1hat), abs(khat&E2hat), abs(E1hat&E2hat), abs(khatlen-1.0), abs(E1len-1.0), abs(E2len-1.0));

  assert( abs(khat&E1hat) < 1E-8 );
  assert( abs(khat&E2hat) < 1E-8 );
  assert( abs(E1hat&E2hat) < 1E-8 );
  assert( abs(khatlen-1.0) < 1E-8 );
  assert( abs(E1len-1.0) < 1E-8 );
  assert( abs(E2len-1.0) < 1E-8 );
}

void CoccolithSimulation::computeEvectorOrthogonalToPropagation( const meep::vec &r, meep::vec &Eperp, meep::vec &Epar )
{
  meep::vec incWave;
  RotationAxis_t rotAx;
  meep::vec EperpInc;
  meep::vec EparInc;
  switch( propagationDir )
  {
    case MainPropDirection_t::X:
      incWave = meep::vec(1.0,0.0,0.0);
      EperpInc   = meep::vec(0.0,1.0,0.0);
      EparInc  = meep::vec(0.0,0.0,1.0);
      rotAx = RotationAxis_t::Z;
      break;
    case MainPropDirection_t::Y:
      incWave = meep::vec(0.0,1.0,0.0);
      EperpInc   = meep::vec(0.0,0.0,1.0);
      EparInc   = meep::vec(1.0,0.0,0.0);
      rotAx = RotationAxis_t::X;
      break;
    case MainPropDirection_t::Z:
      EperpInc   = meep::vec(1.0,0.0,0.0);
      EparInc   = meep::vec(0.0,1.0,0.0);
      incWave = meep::vec(0.0,0.0,1.0);
      rotAx = RotationAxis_t::Y;
      break;
  }

  static const double PI = acos(-1.0);
  if ( srcPos == SourcePosition_t::BOTTOM )
  {
    // Rotate PI around the first axis
    double rotMat[3][3];
    setUpRotationMatrix(rotAx, PI, rotMat );
    incWave = rotateVector( rotMat, incWave );
    EperpInc = rotateVector( rotMat, EperpInc );
    EparInc = rotateVector( rotMat, EparInc );
  }

  Eperp = cross(incWave,r);
  Eperp = Eperp/norm(Eperp);
  Epar = cross(r,Eperp);
  Epar = Epar/norm(Epar);
  assert( abs((Epar&r)/norm(r)) < 1E-6 );
  assert( abs((Eperp&r)/norm(r)) < 1E-6 );
  assert( abs(Eperp&Epar) < 1E-6 );

  // Set the current rotation angle
  currentStokesVectorRotationAngleRad = acos(Eperp&EperpInc);
  if ( (Eperp&EparInc) < 0.0 ) currentStokesVectorRotationAngleRad = -currentStokesVectorRotationAngleRad;
}

void CoccolithSimulation::getLocalStokes( double theta, double phi, const cdouble EH[], LocalStokes &locStoke, Stokes &inc )
{
  locStoke.I.resize(n2fBox->Nfreq);
  locStoke.Q.resize(n2fBox->Nfreq);
  locStoke.U.resize(n2fBox->Nfreq);
  locStoke.V.resize(n2fBox->Nfreq);

  meep::vec E1hat;
  meep::vec E2hat;

  double x = sin(theta)*cos(phi);
  double y = sin(theta)*sin(phi);
  double z = cos(theta);
  permumteToFitPropDir(x,y,z);
  meep::vec rhat(x,y,z);
  computeEvectorOrthogonalToPropagation( rhat, E1hat, E2hat );
  bool EfieldWasStored = false;
  inc = incStoke;
  inc.rotate( currentStokesVectorRotationAngleRad );
  for ( unsigned int i=0;i<locStoke.I.size();i++ )
  {
    cdouble Ex = EH[6*i];
    cdouble Ey = EH[6*i+1];
    cdouble Ez = EH[6*i+2];
    //permumteToFitPropDir(Ex,Ey,Ez);
    cdouble E1 = E1hat.x()*Ex + E1hat.y()*Ey + E1hat.z()*Ez;
    cdouble E2 = E2hat.x()*Ex + E2hat.y()*Ey + E2hat.z()*Ez;
    locStoke.I[i] = pow(abs(E1),2) + pow(abs(E2),2);
    locStoke.Q[i] = pow(abs(E1),2) - pow(abs(E2),2);
    locStoke.U[i] = 2.0*(E1*conj(E2)).real();
    locStoke.V[i] = -2.0*(E1*conj(E2)).imag();
    //locStoke.V[i] = 2.0*(E1*conj(E2)).imag();
    // Debug: Make sure that there is no field component in the radial direction
    double intensityPhiTheta = pow(abs(E1),2) + pow(abs(E2),2);
    double intensityXYZ = pow(abs(EH[6*i]),2) + pow(abs(EH[6*i+1]),2) + pow(abs(EH[6*i+2]),2);
    double relDev = abs(intensityPhiTheta - intensityXYZ )/intensityXYZ;
    if ( relDev > relativeRadFieldCompThreshold ) numberOfTimesThereIsRadialFieldComponent++;
    totalNumberOfFarFieldEvaluations++;

    if ( i==n2fBox->Nfreq/2 )
    {
      // Store the electric field for the center frequency only
      locStoke.Ephi = abs(E2);
      locStoke.Etheta = abs(E1);
      EfieldWasStored = true;
    }
  }
  assert( EfieldWasStored );
}

void CoccolithSimulation::loadBoundingCurrents( const char* fname )
{
  if ( field == NULL )
  {
    throw (runtime_error("Field needs to be initialized before loading currents!") );
  }
  else if ( n2fBox == NULL )
  {
    throw ( runtime_error("Near to far field box needs to be initialized!"));
  }
  n2fBox->load_hdf5( *field, fname );
}

void CoccolithSimulation::saveStokesPhiTheta()
{
  assert( file != NULL );
  for ( unsigned int i=0;i<stokesI.size();i++ )
  {
    int nrow = stokesI[i].n_rows;
    int ncol = stokesI[i].n_cols;
    int rank[2] = {nrow, ncol};
    arma::mat stokes = stokesI[i].t();
    {
      stringstream dset;
      dset << "stokesI" << i;
      file->write( dset.str().c_str(), 2, rank, stokes.memptr(), false );
      file->prevent_deadlock();
    }

    {
      stringstream dset;
      dset << "stokesQ" << i;
      stokes = stokesQ[i].t();
      file->write( dset.str().c_str(), 2, rank, stokes.memptr(), false );
      file->prevent_deadlock();
    }

    {
      stringstream dset;
      dset << "stokesU" << i;
      stokes = stokesU[i].t();
      file->write( dset.str().c_str(), 2, rank, stokes.memptr(), false );
      file->prevent_deadlock();
    }

    {
      stringstream dset;
      dset << "stokesV" << i;
      stokes = stokesV[i].t();
      file->write( dset.str().c_str(), 2, rank, stokes.memptr(), false );
      file->prevent_deadlock();
    }
  }

  int nrows = Ephi.n_rows;
  int ncols = Ephi.n_cols;
  int rank[2] = {nrows,ncols};
  arma::mat Etrans = Ephi.t();
  file->write("Ephi", 2, rank, Etrans.memptr(), false );
  file->prevent_deadlock();

  Etrans = Etheta.t();
  file->write("Etheta", 2, rank, Etrans.memptr(), false );
  file->prevent_deadlock();

  int size = stokesIInc.size();
  file->write("StokesIInc", 1, &size, &stokesIInc[0], false );
  file->prevent_deadlock();

  file->write("StokesQInc", 1, &size, &stokesQInc[0], false );
  file->prevent_deadlock();

  file->write("StokesUInc", 1, &size, &stokesUInc[0], false );
  file->prevent_deadlock();

  file->write("StokesVInc", 1, &size, &stokesVInc[0], false );
  file->prevent_deadlock();

  meep::master_printf("In %d of %d of the field evaluations there seem to be a radial field component", numberOfTimesThereIsRadialFieldComponent, totalNumberOfFarFieldEvaluations );
}

void CoccolithSimulation::setIncStokesVector( const int stk[4] )
{
  incStoke.I = stk[0];
  incStoke.Q = stk[1];
  incStoke.U = stk[2];
  incStoke.V = stk[3];
  bool isASupportedStokesVector = false;
  for (unsigned int i=0;i<supportedStokes.size();i++ )
  {
    if ( incStoke == supportedStokes[i])
    {
      isASupportedStokesVector = true;
      break;
    }
  }

  if ( !isASupportedStokesVector )
  {
    meep::abort("The supplied stokes vector is not supported!");
  }
  initialStokesVectorSet = true;
}

meep::vec CoccolithSimulation::cross( const meep::vec &v1, const meep::vec &v2 )
{
  double x = v1.y()*v2.z() - v1.z()*v2.y();
  double y = v1.z()*v2.x() - v1.x()*v2.z();
  double z = v1.x()*v2.y() - v1.y()*v2.x();
  return meep::vec(x,y,z);
}

double CoccolithSimulation::norm( const meep::vec &vec )
{
  return sqrt( pow(vec.x(),2) + pow(vec.y(),2) + pow(vec.z(),2) );
}

void CoccolithSimulation::setUpRotationMatrix( RotationAxis_t raxis, double alpha, double matrix[3][3] )
{

  switch(raxis)
  {
    case RotationAxis_t::X:
      matrix[0][0] = 1.0;
      matrix[0][1] = 0.0;
      matrix[0][2] = 0.0;
      matrix[1][0] = 0.0;
      matrix[1][1] = cos(alpha);
      matrix[1][2] = sin(alpha);
      matrix[2][0] = 0.0;
      matrix[2][1] = -sin(alpha);
      matrix[2][2] = cos(alpha);
      break;
    case RotationAxis_t::Y:
      matrix[0][0] = cos(alpha);
      matrix[0][1] = 0.0;
      matrix[0][2] = sin(alpha);
      matrix[1][0] = 0.0;
      matrix[1][1] = 1.0;
      matrix[1][2] = 0.0;
      matrix[2][0] = -sin(alpha);
      matrix[2][1] = 0.0;
      matrix[2][2] = cos(alpha);
      break;
    case RotationAxis_t::Z:
      matrix[0][0] = cos(alpha);
      matrix[0][1] = sin(alpha);
      matrix[0][2] = 0.0;
      matrix[1][0] = -sin(alpha);
      matrix[1][1] = cos(alpha);
      matrix[1][2] = 0.0;
      matrix[2][0] = 0.0;
      matrix[2][1] = 0.0;
      matrix[2][2] = 1.0;
      break;
  }
}

meep::vec CoccolithSimulation::rotateVector( const double rotMat[3][3], const meep::vec &vec )
{
  double x = rotMat[0][0]*vec.x() + rotMat[0][1]*vec.y() + rotMat[0][2]*vec.z();
  double y = rotMat[1][0]*vec.x() + rotMat[1][1]*vec.y() + rotMat[1][2]*vec.z();
  double z = rotMat[2][0]*vec.x() + rotMat[2][1]*vec.y() + rotMat[2][2]*vec.z();
  return meep::vec(x,y,z);
}

void CoccolithSimulation::combineRotationMatrices( const double firstRot[3][3], const double secondRot[3][3], double combined[3][3] )
{
  for ( unsigned int i=0;i<3;i++ )
  for ( unsigned int j=0;j<3;j++ )
  for ( unsigned int k=0;k<3;k++ )
  {
    combined[i][j] += secondRot[i][k]*firstRot[k][j];
  }
}

void CoccolithSimulation::printRotationMatrix( const double rotMat[3][3] )
{
  cout <<"--------------------------\n";
  for ( unsigned int i=0;i<3;i++ )
  {
    for ( unsigned int j=0;j<3;j++ )
    {
      //if ( meep::am_master() )
      {
        cout << rotMat[i][j] << " ";
      }
    }
    cout << endl;
  }
  cout <<"--------------------------\n";
}

void CoccolithSimulation::addIdentifierToBackups( const char* extra )
{
  n2fBoxBackup += extra;
  reflFluxBoxBackup += extra;
  reflFluxPlaneBackup += extra;
}
//============================ STOKES CLASS ====================================

bool Stokes::operator==( const Stokes &other ) const
{
  double EPS=1E-10;
  return (abs(this->I-other.I)<EPS) && (abs(this->Q-other.Q)<EPS) && \
  (abs(this->U-other.U)<EPS) && (abs(this->V-other.V)<EPS);
}

void Stokes::rotate( double angleRad )
{
  double Q0 = Q;
  double U0 = U;
  Q = cos(2.0*angleRad)*Q0 + sin(2.0*angleRad)*U0;
  U = cos(2.0*angleRad)*U0 - sin(2.0*angleRad)*Q0;
}
