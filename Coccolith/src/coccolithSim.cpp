#include "coccolithSim.hpp"
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <omp.h>
#include <ctime>
#include <iomanip>

using namespace std;

meep::vec CoccolithSimulation::waveVec;

// Only for debugging
double sigmaTest( const meep::vec &r )
{
  return 1.0;
}

CoccolithSimulation::~CoccolithSimulation()
{
  if ( srcVol != NULL ) delete srcVol; srcVol=NULL;
  // Sources and flux are deleted in the destructor of field

  // Apperently it looks like field needs to be deleted before struc. Bug??
  if ( field != NULL ) delete field; field=NULL;
  if ( struc != NULL ) delete struc; struc=NULL;
  if ( dftVolTransmit != NULL ) delete dftVolTransmit; dftVolTransmit=NULL;
  delete dftVolRefl; dftVolRefl = NULL;
  delete dftVolBox; dftVolBox = NULL;
  delete transmitFlux; transmitFlux = NULL;
  delete reflFlux; reflFlux = NULL;
  delete fluxBox; fluxBox = NULL;
  delete file; file=NULL;

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
  meep::vec reflDFTCrn1;
  meep::vec reflDFTCrn2;

  // Coordinates for flux box
  meep::vec boxCrn1 = meep::vec( gdvol.xmin()+getPMLThickness(), gdvol.ymin()+getPMLThickness(), gdvol.zmin()+getPMLThickness() );
  meep::vec boxCrn2 = meep::vec( gdvol.xmax()-getPMLThickness(), gdvol.ymax()-getPMLThickness(), gdvol.zmax()-getPMLThickness() );

  switch ( propagationDir )
  {
    case MainPropDirection_t::X:
      transDFTCrn1 = meep::vec( getTransFluxPos(), gdvol.ymin(), gdvol.zmin() );
      transDFTCrn2 = meep::vec( getTransFluxPos(), gdvol.ymax(), gdvol.zmax() );
      reflDFTCrn1 = meep::vec( getSrcFluxPos(), gdvol.ymin(), gdvol.zmin() );
      reflDFTCrn2 = meep::vec( getSrcFluxPos(), gdvol.ymax(), gdvol.zmax() );
      break;
    case MainPropDirection_t::Y:
      transDFTCrn1 = meep::vec( gdvol.xmin(), getTransFluxPos(), gdvol.zmin() );
      transDFTCrn2 = meep::vec( gdvol.xmax(), getTransFluxPos(), gdvol.zmax() );
      reflDFTCrn1 = meep::vec( gdvol.xmin(), getSrcFluxPos(), gdvol.zmin() );
      reflDFTCrn2 = meep::vec( gdvol.xmax(), getSrcFluxPos(), gdvol.zmax() );
      break;
    case MainPropDirection_t::Z:
      transDFTCrn1 = meep::vec( gdvol.xmin(), gdvol.ymin(), getTransFluxPos() );
      transDFTCrn2 = meep::vec( gdvol.xmax(), gdvol.ymax(), getTransFluxPos() );
      reflDFTCrn1 = meep::vec( gdvol.xmin(), gdvol.ymin(), getSrcFluxPos() );
      reflDFTCrn2 = meep::vec( gdvol.xmax(), gdvol.ymax(), getSrcFluxPos() );
      break;
  }

  if ( dftVolTransmit != NULL ) delete dftVolTransmit;
  dftVolTransmit = new meep::volume( transDFTCrn1, transDFTCrn2 );

  transmitFlux = new meep::dft_flux( field->add_dft_flux_plane( *dftVolTransmit, centerFrequency-freqWidth/2.0, centerFrequency+freqWidth/2.0, nfreq ) );

  dftVolRefl = new meep::volume( reflDFTCrn1, reflDFTCrn2 );
  dftVolBox = new meep::volume( boxCrn1, boxCrn2 );
  reflFlux = new meep::dft_flux( field->add_dft_flux_plane( *dftVolRefl, centerFrequency-freqWidth/2.0, centerFrequency+freqWidth/2.0, nfreq ) );
  fluxBox = new meep::dft_flux( field->add_dft_flux_box( *dftVolBox, centerFrequency-freqWidth/2.0, centerFrequency+freqWidth/2.0, nfreq ) );

  if ( !material->isReferenceRun() )
  {
    // There should be a backup file with the fluxes stored from the reference run
    reflFlux->load_hdf5( *field, reflFluxPlaneBackup.c_str() );
    fluxBox->load_hdf5( *field, reflFluxBoxBackup.c_str() );
    reflFlux->scale_dfts(-1.0);
    fluxBox->scale_dfts(-1.0);
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
  ss << "voxelMaterialSim_" << uid; // File extension added automatically
  if ( file == NULL )
  {
    if ( material->isReferenceRun() )
    {
      file = field->open_h5file( ss.str().c_str() );
      saveGeometry();
      saveDFTParameters();
    }
    else
    {
      file = field->open_h5file( ss.str().c_str(), meep::h5file::READWRITE );
    }
  }
  if ( !material->isReferenceRun() )
  {
    field->output_hdf5( meep::Dielectric, gdvol.surroundings(), file, false, true );
    file->prevent_deadlock();
  }
  saveDFTSpectrum();
  //saveDFTStokes();

  if (( meep::am_master() ) && (!material->isReferenceRun() ))
  {
    clog << "Results written to " << ss.str() << endl;
  }

  /*
  if (( meep::am_master() ) && (!material->isReferenceRun() ))
  {
      clog << "Saving parameters...\n";
      saveDFTParameters();
      saveGeometry();
      Json::StyledWriter sw;
      ofstream of;
      string fname = outdir+"/"+prefix+"-parameters.json";
      of.open( fname.c_str() );
      if ( !of.good() )
      {
        clog << "Could not open JSON file " << fname << endl;
        return;
      }
      of << sw.write(root) << endl;
      of.close();
  */
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
  }
  else
  {
    int length = transmitFlux->Nfreq;
    file->write( "spectrumTransmitted", 1, &length, transmitFlux->flux(), false );
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

  /*
  Json::Value dftParams;
  dftParams["freqmin"] = transmitFlux->freq_min;
  dftParams["dfreq"] = transmitFlux->dfreq;
  dftParams["Nfreq"] = transmitFlux->Nfreq;

  // Write geometrical positions
  meep::vec mincrn = dftVolTransmit->get_min_corner();
  meep::vec maxcrn = dftVolTransmit->get_max_corner();
  Json::Value corner(Json::arrayValue);
  corner.append( mincrn.x() );
  corner.append( mincrn.y() );
  corner.append( mincrn.z() );
  dftParams["corner1"] = corner;

  Json::Value corner2(Json::arrayValue);
  corner2.append( maxcrn.x() );
  corner2.append( maxcrn.y() );
  corner2.append( maxcrn.z() );
  dftParams["corner2"] = corner2;
  root["dft"] = dftParams;
  */
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

  /*
  //assert( file != NULL );
  Json::Value geometry;

  meep::vec mincrn = srcVol->get_min_corner();
  meep::vec maxcrn = srcVol->get_max_corner();

  Json::Value corner1(Json::arrayValue);
  corner1.append( mincrn.x() );
  corner1.append( mincrn.y() );
  corner1.append( mincrn.z() );
  geometry["corner1"] = corner1;
  Json::Value corner2(Json::arrayValue);
  corner2.append( maxcrn.x() );
  corner2.append( maxcrn.y() );
  corner2.append( maxcrn.z() );
  geometry["corner2"] = corner2;

  Json::Value waveVector(Json::arrayValue);
  waveVector.append( waveVec.x() );
  waveVector.append( waveVec.y() );
  waveVector.append( waveVec.z() );
  geometry["wavVector"] = waveVector;

  int single = 1;
  geometry["wavelength"] = getWavelength();

  int nFreqParam = 2;
  double freqParam[2] = {centerFrequency, freqWidth};
  geometry["dimLessAngularFreq"] = centerFrequency;
  geometry["dimLessFreqWidth"] = freqWidth;
  geometry["xmin"] = gdvol.xmin();
  geometry["xmax"] = gdvol.xmax();
  geometry["ymin"] = gdvol.ymin();
  geometry["ymax"] = gdvol.ymax();
  geometry["zmin"] = gdvol.zmin();
  geometry["zmax"] = gdvol.zmax();

  geometry["pmlThickness"] = getPMLThickness();

  geometry["voxelSize"] = material->getVoxelSize();
  root["geometry"] = geometry;*/
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

void CoccolithSimulation::saveDFTStokes()
{
  if ( transmitFlux == NULL )
  {
    throw( runtime_error("Flux plane needs to be set before saving DFT stokes parameters!") );
  }
  stokes.compute( transmitFlux->E );
  int length = stokes.getI().size();
  if ( material->isReferenceRun() )
  {
    file->write( "stokeRefI", 1, &length, &stokes.getI()[0], false );
    file->prevent_deadlock();
    file->write( "stokeRefQ", 1, &length, &stokes.getQ()[0], false );
    file->prevent_deadlock();
    file->write( "stokeRefU", 1, &length, &stokes.getU()[0], false );
    file->prevent_deadlock();
    file->write( "stokeRefV", 1, &length, &stokes.getV()[0], false );
    file->prevent_deadlock();
  }
  else
  {
    file->write( "stokeTransI", 1, &length, &stokes.getI()[0], false );
    file->prevent_deadlock();
    file->write( "stokeTransQ", 1, &length, &stokes.getQ()[0], false );
    file->prevent_deadlock();
    file->write( "stokeTransU", 1, &length, &stokes.getU()[0], false );
    file->prevent_deadlock();
    file->write( "stokeTransV", 1, &length, &stokes.getV()[0], false );
    file->prevent_deadlock();
  }

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
