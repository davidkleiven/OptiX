#include "genericScattering.hpp"
#include <stdexcept>

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
  ef.setExportDimensions( exportNx, exportNy );
  ei.setExportDimensions( exportNx, exportNy );
  ep.setExportDimensions( exportNx, exportNy );
  ff.setPadLength( FFTPadLength );
  ff.setReference( reference );

  setTransverseDiscretization( xmin, xmax, dx, downSampleX );
  setVerticalDiscretization( ymin, ymax, dy );
  setLongitudinalDiscretization( zmin, zmax, (zmax-zmin)/3.0, 1 ); // Settings for reference run
  setWaveLength( wavelength );
  gbeam.setWavelength( wavelength );
  setBoundaryConditions( gbeam );
  fft3Dsolver.setIntensityMinMax( intensityMin, intensityMax );
  fft3Dsolver.setPhaseMinMax( phaseMin, phaseMax );

  if ( imgname != "" )
  {
    fft3Dsolver.storeImages( imgname.c_str() );
  }

  setSolver( fft3Dsolver );
  *this << ef << ei << ep << ff;
}

void GenericScattering::getXrayMatProp( double x, double y, double z, double &delta, double &beta )
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
  clog << "Reference solution computed\n";

  reset();

  fft3Dsolver.visualizeRealSpace();
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
