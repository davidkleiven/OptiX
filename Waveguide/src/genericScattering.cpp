#include "genericScattering.hpp"
#include <stdexcept>
//#define PRINT_DEBUG

using namespace std;
GenericScattering::GenericScattering( const char* name ): ParaxialSimulation(name)
{
  gbeam.setDim( ParaxialSource::Dim_t::THREE_D );
}

void GenericScattering::setMaxScatteringAngle( double anglemax )
{
  ff.setAngleRange( -anglemax, anglemax );
}

void GenericScattering::init()
{
  #ifdef PRINT_DEBUG
    clog << "Initialization started...\n";
    clog << "Set export dimensions and padlength...\n";
  #endif

  ef.setExportDimensions( exportNx, exportNy );
  ei.setExportDimensions( exportNx, exportNy );
  ep.setExportDimensions( exportNx, exportNy );
  ff.setPadLength( FFTPadLength );

  #ifdef PRINT_DEBUG
    clog << "Set reference solution array...\n";
  #endif

  #ifdef PRINT_DEBUG
    clog << "Set discretization...\n";
  #endif

  setTransverseDiscretization( xmin, xmax, dx, downSampleX );
  setVerticalDiscretization( ymin, ymax, dy );
  setLongitudinalDiscretization( zmin, zmax, (zmax-zmin)/3.0, 1 ); // Settings for reference run

  #ifdef PRINT_DEBUG
    clog << "Set wavelength...\n";
  #endif

  setWaveLength( wavelength );
  gbeam.setWavelength( wavelength );

  #ifdef PRINT_DEBUG
    clog << "Set visualization phase min/max and intensity min/max...\n";
  #endif

  if ( realTimeVisualization )
  {
    fft3Dsolver.visualizeRealSpace();
    fft3Dsolver.setIntensityMinMax( intensityMin, intensityMax );
    fft3Dsolver.setPhaseMinMax( phaseMin, phaseMax );
  }

  if ( imgname != "" )
  {
    fft3Dsolver.storeImages( imgname.c_str() );
  }

  #ifdef PRINT_DEBUG
    clog << "Set solver and post processing modules...\n";
  #endif

  setSolver( fft3Dsolver );

  #ifdef PRINT_DEBUG
    clog << "Set boundary conditions from Gaussian beam...\n";
  #endif

  setBoundaryConditions( gbeam );

  *this << ef << ei << ep << ff;

  #ifdef PRINT_DEBUG
    clog << "Initialization finished...\n";
  #endif
}

void GenericScattering::getXrayMatProp( double x, double y, double z, double &delta, double &beta ) const
{
  if ( isReferenceRun )
  {
    delta = 0.0;
    beta = 0.0;
    return;
  }
  material->getXrayMatProp( x, y, z, delta, beta );
}

void GenericScattering::solve()
{
  if ( material == NULL )
  {
    throw( runtime_error("No material set!") );
  }
  init();
  printInfo();

  // Reference run
  ParaxialSimulation::solve();

  reference = fft3Dsolver.getLastSolution3D();
  ff.setReference( reference );
  clog << "Reference solution computed\n";

  reset();
  isReferenceRun = false;

  if ( realTimeVisualization )
  {
    fft3Dsolver.visualizeRealSpace();
  }
  // Set resolution for higher
  setLongitudinalDiscretization( zmin, zmax, dz, downSampleZ );
  ParaxialSimulation::solve();
}

void GenericScattering::printInfo() const
{
  cout << "===================== SIMULATION INFO ===========================\n";
  cout << "Maximum size of matrices for export: Nx= " << exportNx << " Ny=" << exportNy << endl;
  cout << "FFT pad length: " << FFTPadLength << endl;
  cout << "Domain size: x=[" << xmin << "," << xmax << "], y=[" << ymin << "," << ymax << "], z=["<<zmin<<","<<zmax<<"]\n";
  cout << "Stepsize: dx=" << dx << ", dy="<<dy<<", dz="<<dz<<endl;
  cout << "Downsample 3D solution by factor: downX=" << downSampleX << " downY=" << downSampleY << " downZ="<<downSampleZ << endl;
  cout << "Wavelength: " << wavelength << endl;
  cout << "Intensity limit for realtime visualization: Imin=" << intensityMin << " Imax=" << intensityMax << endl;
  cout << "Phase limit for realtime visualization: Pmin=" << phaseMin << " Pmax=" << phaseMax << endl;
  cout << "Image are saved to: " << imgname << endl;
  cout << "Gaussian beam waist: " << gbeam.getWaist() << endl;
  cout << "Scattering angle: [" << ff.getMinScatteringAngle() << "," << ff.getMaxScatteringAngle() << "]\n";
  cout << "==================================================================\n";
}
