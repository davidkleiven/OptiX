#include "coccolithSim.hpp"
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <omp.h>
#include <ctime>

using namespace std;

meep::vec CoccolithSimulation::waveVec;

CoccolithSimulation::~CoccolithSimulation()
{
  // Sources and flux are deleted in the destructor of field
  if ( srcVol != NULL ) delete srcVol; srcVol=NULL;
  if ( struc != NULL ) delete struc; struc=NULL;
  if ( field != NULL ) delete field; field=NULL;
  if ( dftVolTransmit != NULL ) delete dftVolTransmit; dftVolTransmit=NULL;
  if ( source != NULL ) delete source; source=NULL;

  if ( monitor1 != NULL ) delete monitor1; monitor1=NULL;
  if ( monitor2 != NULL ) delete monitor2; monitor2=NULL;
  #ifdef HAVE_LIB_VISA
    if ( plots != NULL ) delete plots; plots=NULL;
  #endif
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
      break;
    case MainPropDirection_t::Y:
      fieldComp = meep::component::Ez;
      break;
    case MainPropDirection_t::Z:
      fieldComp = meep::component::Ex;
      break;
  }

  field->add_volume_source( fieldComp, *sourceTime, *srcVol, amplitude );
}

void CoccolithSimulation::init()
{
  if ( material == NULL )
  {
    throw( runtime_error("No material set! Needs to be set before call init()!") );
  }
  domainInfo();
  initializeGeometry();

  material->setDomainSize( gdvol, getPMLThickness() );
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

  switch ( propagationDir )
  {
    case MainPropDirection_t::X:
      transDFTCrn1 = meep::vec( getTransFluxPos(), gdvol.ymin(), gdvol.zmin() );
      transDFTCrn2 = meep::vec( getTransFluxPos(), gdvol.ymax(), gdvol.zmax() );
      break;
    case MainPropDirection_t::Y:
      transDFTCrn1 = meep::vec( gdvol.xmin(), getTransFluxPos(), gdvol.zmin() );
      transDFTCrn2 = meep::vec( gdvol.xmax(), getTransFluxPos(), gdvol.zmax() );
      break;
    case MainPropDirection_t::Z:
      transDFTCrn1 = meep::vec( gdvol.xmin(), gdvol.ymin(), getTransFluxPos() );
      transDFTCrn2 = meep::vec( gdvol.xmax(), gdvol.ymax(), getTransFluxPos() );
      break;
  }

  if ( dftVolTransmit != NULL ) delete dftVolTransmit;
  dftVolTransmit = new meep::volume( transDFTCrn1, transDFTCrn2 );

  if ( transmitFlux != NULL ) delete transmitFlux;
  transmitFlux = new meep::dft_flux( field->add_dft_flux_plane( *dftVolTransmit, centerFrequency-freqWidth/2.0, centerFrequency+freqWidth/2.0, nfreq ) );
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
      return getSrcPos() + 5.0*material->getVoxelSize();
    case SourcePosition_t::BOTTOM:
      return getSrcPos() - 5.0*material->getVoxelSize();
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

  if ( meep::my_rank() == 0 )
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

  if ( meep::my_rank() == 0 )
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

  if ( meep::my_rank() == 0 )
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
  if ( meep::my_rank() == 0 )
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
    stringstream ss;
    ss << localtime(&t)->tm_year+1900 << localtime(&t)->tm_mon << localtime(&t)->tm_mday <<"_"<<
    localtime(&t)->tm_hour <<localtime(&t)->tm_min << localtime(&t)->tm_sec;
    uid = ss.str();
  }
}
void CoccolithSimulation::exportResults()
{
  setUID();

  if ( material == NULL )
  {
    throw( runtime_error("Material needs to be set before calling exportResults()!") );
  }

  stringstream ss;
  ss << "data/voxelMaterialSim_" << uid; // File extension added automatically
  if ( file == NULL )
  {
    file = field->open_h5file( ss.str().c_str() );
    saveGeometry();
  }

  saveDFTSpectrum();

  if ( !material->isReferenceRun() )
  {
    field->output_hdf5( meep::Dielectric, gdvol.surroundings(), file, false, true );
  }
  if ( meep::my_rank() == 0 )
  {
      clog << "Results written to " << ss.str() << endl;
  }
}

void CoccolithSimulation::saveDFTSpectrum()
{
  assert( file != NULL );
  assert( material != NULL );

  if ( material->isReferenceRun() )
  {
    int length = transmitFlux->Nfreq;
    file->write( "spectrumReference", 1, &length, transmitFlux->flux(), false );
    return;
  }
  else
  {
    int length = transmitFlux->Nfreq;
    file->write( "spectrumTransmitted", 1, &length, transmitFlux->flux(), false );
  }

  const int nParams = 3;
  double stat[nParams];
  stat[0] = transmitFlux->freq_min;
  stat[1] = transmitFlux->dfreq;
  stat[2] = transmitFlux->Nfreq;

  // Write parameters
  file->write( "spectrumFreqs", 1, &nParams, stat, false );

  // Write geometrical positions
  meep::vec mincrn = dftVolTransmit->get_min_corner();
  meep::vec maxcrn = dftVolTransmit->get_max_corner();
  int nCrns = 6;
  double cornersTrans[nCrns] = {mincrn.x(), mincrn.y(), mincrn.z(), maxcrn.x(), maxcrn.y(), maxcrn.z()};
  file->write( "spectrumCornersTransmit", 1, &nCrns, cornersTrans, false );
}

void CoccolithSimulation::saveGeometry()
{
  assert( file != NULL );

  meep::vec mincrn = srcVol->get_min_corner();
  meep::vec maxcrn = srcVol->get_max_corner();
  int nCrn = 6;
  double corners[6] = {mincrn.x(), mincrn.y(), mincrn.z(), maxcrn.x(), maxcrn.y(), maxcrn.z()};
  file->write( "geometrySourceVolume", 1, &nCrn, corners, false );

  int ndim = 3;
  double incWaveVec[ndim] = {waveVec.x(), waveVec.y(), waveVec.z()};
  file->write( "geometryIncWavevector", 1, &ndim, incWaveVec, false );

  int single = 1;
  double wavelength = getWavelength();
  file->write( "geometryWavelength", 1, &single, &wavelength, false );

  int nFreqParam = 2;
  double freqParam[2] = {centerFrequency, freqWidth};
  file->write( "geometryDimlessAngularFreq", 1, &nFreqParam, freqParam, false );

  double domainSize[nCrn] = {gdvol.xmin(), gdvol.xmax(), gdvol.ymin(), gdvol.ymax(), gdvol.zmin(), gdvol.zmax()};
  file->write( "geometryDomainSize", 1, &nCrn, domainSize, false );

  double pmlThickness = getPMLThickness();
  file->write( "geometryPmlThickness", 1, &single, &pmlThickness, false );

  double vxsize = material->getVoxelSize();
  file->write( "voxelsizeInNanoMeter", 1, &single, &vxsize, false );
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
  gdvol = meep::vol3d( material->sizeX()+2.0*pmlThickPx, material->sizeY()+2.0*pmlThickPx, material->sizeZ()+2.0*pmlThickPx, resolution );
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
