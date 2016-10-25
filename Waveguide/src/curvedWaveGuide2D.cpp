#include "curvedWaveGuide2D.hpp"
#include "waveGuide.hpp"
#include "solver1D.hpp"
#include <cassert>
#include <cmath>
#include <vector>
#include "cladding.hpp"
#include "solver2D.hpp"
#include <cmath>
#include "controlFile.hpp"
#include "linearMap1D.hpp"
#include <H5Cpp.h>
#include <hdf5_hl.h>
//#include "source.hpp"
#include "gaussianBeam.hpp"
#include "planeWave.hpp"

#define PROJECTION_DEBUG

using namespace std;

bool CurvedWaveGuideFD::isInsideGuide( double x, double z ) const
{
//  double d = sqrt(x*x + z*z);
  return (2.0*x*R+z*z > 0.0 ) && (2.0*x*R+z*z < 2.0*width*R);
  //return (pow(R+x,2)+z*z > pow(R,2)) && (pow(R+x,2)+z*z < pow(R+width,2));
  //return ( d < R+width) && (d > R );
}

void CurvedWaveGuideFD::fillInfo( Json::Value &obj ) const
{
  obj["RadiusOfCurvature"] = R;
  obj["Width"] = width;
  obj["crd"] = "cartesian";
}

void CurvedWaveGuideFD::computeTransmission( double step )
{
  // Integrate across waveguide
  double fluxAtZero = 0.0;
  unsigned int wgStart = 0;
  unsigned int wgEnd = 0;
  unsigned int zIndx;
  closestIndex( 0.0, 0.0, wgStart, zIndx );
  closestIndex( width, 0.0, wgEnd, zIndx );
  double intensityAtZero = trapezoidalIntegrateIntensityZ( zIndx, wgStart, wgEnd );

  // Just simple trapezoidal rule for integrating in z-direction
  //double z = step;
  //while ( z < zDisc->max-step )
  for ( unsigned int iz=0;iz<nodeNumberLongitudinal();iz++ )
  {
    double z = zDisc->min +iz*zDisc->step;
    double xWgStart = waveGuideStartX( z );
    double xWgEnd = waveGuideEndX( z );
    closestIndex( xWgStart, z, wgStart, zIndx );
    closestIndex( xWgEnd, z, wgEnd, zIndx );

    // Assertions for debugging (allow first and last point to be outside )
    assert( isInsideGuide( xWgStart+1, z ) );
    assert( isInsideGuide( xWgEnd-1, z) );

    double intensity = trapezoidalIntegrateIntensityZ( zIndx, wgStart, wgEnd );
    transmission.push_back( intensity/intensityAtZero );
    z += step;
  }
  stepWhenComputingTransmission = step; // Save for later
}

double CurvedWaveGuideFD::waveGuideStartX( double z ) const
{
  return -0.5*z*z/R;
}

double CurvedWaveGuideFD::waveGuideEndX( double z ) const
{
  return waveGuideStartX( z ) + width;
}

void CurvedWaveGuideFD::saveTransmission( ControlFile &ctl ) const
{
  Json::Value trans;
  trans["zStart"] = zDisc->step;
  trans["zEnd"] = zDisc->max;
  trans["step"] = stepWhenComputingTransmission;
  string fname = ctl.getFnameTemplate();
  fname += "_trans.h5";
  trans["file"] = fname;

  int rank = 1;
  hsize_t dim = transmission.size();
  hid_t file_id = H5Fcreate( fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  H5LTmake_dataset( file_id, "transmission", rank, &dim, H5T_NATIVE_DOUBLE, &transmission[0]);
  H5Fclose(file_id);
  ctl.get()["Transmission"] = trans;
  clog << "Transmission is written to " << fname << endl;
}

void CurvedWaveGuideFD::init( const ControlFile &ctl )
{
  WaveGuideFDSimulation::init(ctl);
  width = ctl.get()["waveguide"]["Width"].asDouble();
  R = ctl.get()["waveguide"]["RadiusOfCurvature"].asDouble();
}

void CurvedWaveGuideFD::getFieldInsideWG( arma::mat &matrix ) const
{
  unsigned int nNodesTransverseWG = width/xDisc->step;

  matrix.set_size( nNodesTransverseWG, nodeNumberLongitudinal() );
  for ( unsigned int iz=0;iz<nodeNumberLongitudinal();iz++ )
  {
    double z = zDisc->min + zDisc->step;
    double x0 = waveGuideStartX( z );
    unsigned int ixStart, dummy;
    closestIndex( x0, z, ixStart, dummy );
    for ( unsigned int ix=0;ix<nNodesTransverseWG; ix++ )
    {
      x0 += xDisc->step;
      unsigned int currentIx;
      assert ( isInsideGuide(x0, z) );
      closestIndex( x0, z, currentIx, dummy );
      matrix(ix, iz) = solver->getSolution()(currentIx, iz).real();
    }
  }
}

double CurvedWaveGuideFD::project( double z, const WaveGuide1DSimulation &eig, unsigned int eigenmode ) const
{
  // Use simple trapezoidal integration scheme with the accuracy of the FD simulation
  double outsideWGReg = 0.2*width;
  double xStart = waveGuideStartX( z );
  double xEnd = waveGuideEndX( z );

  // Define a linear map between the two coordinate systems
  /**
  * The conformal tranformation used maps the inner edge onto x = -width
  * and the outer edge onto x = 0.0
  */
  LinearMap1D map;
  CorrespondingPoints p1;
  CorrespondingPoints p2;
  p1.xFrom = xStart;
  p1.xTo = -eig.getWidth();
  p2.xFrom = xEnd;
  p2.xTo = 0.0;
  map.initialize( p1, p2 );

  // Add some buffer to outside the wg
  xStart -= outsideWGReg;
  xEnd += outsideWGReg;

  double dx = xDisc->step;
  unsigned int N = (xEnd-xStart)/dx-1;
  unsigned int xStartIndx, xEndIndx, iz;
  closestIndex( xStart, z, xStartIndx, iz );
  closestIndex( xEnd, z, xEndIndx, iz);
  double transXStart = map.get( xStart );
  double transXend = map.get( xEnd );
  double fval = eig.getSolver()->getSolution(transXStart, eigenmode)*solver->getSolution()(xStartIndx, iz).real();
  double integral = fval;
  fval = eig.getSolver()->getSolution(transXend, eigenmode)*solver->getSolution()(xEndIndx, iz).real();
  integral += fval;

  #ifdef PROJECTION_DEBUG
    vector<double> eigfunc;
    vector<double> fdsolution;
    vector<double> product;
    vector<double> transformedXVec;
  #endif

  for ( unsigned int ix=1;ix<N-1; ix++ )
  {
    double x = xStart + dx*ix;
    double transformedX = map.get( x );
    //cerr << x << " --> " << transformedX << " " <<  eig.getSolver()->getSolution(transformedX, eigenmode) << endl;
    fval = 2.0*eig.getSolver()->getSolution(transformedX, eigenmode)*solver->getSolution()(ix,iz).real();
    integral += fval;

    #ifdef PROJECTION_DEBUG
      eigfunc.push_back( eig.getSolver()->getSolution(transformedX, eigenmode) );
      fdsolution.push_back( solver->getSolution()(ix,iz).real() );
      product.push_back( eig.getSolver()->getSolution(transformedX, eigenmode)*solver->getSolution()(ix,iz).real() );
      transformedXVec.push_back( transformedX );
    #endif
  }

  #ifdef PROJECTION_DEBUG
    int izToDump = 1500;
    int modeToDump = 1;
    if (( iz == izToDump ) && ( eigenmode == modeToDump ) )
    {
      string fname ("data/projectionDebug.h5");
      hid_t file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
      hsize_t dim = eigfunc.size();
      H5LTmake_dataset( file_id, "eigenfunc", 1, &dim, H5T_NATIVE_DOUBLE, &eigfunc[0] );
      H5LTmake_dataset( file_id, "fdsolution", 1, &dim, H5T_NATIVE_DOUBLE, &fdsolution[0] );
      H5LTmake_dataset( file_id, "product", 1, &dim, H5T_NATIVE_DOUBLE, &product[0] );
      H5LTmake_dataset( file_id, "transformedX", 1, &dim, H5T_NATIVE_DOUBLE, &transformedXVec[0] );

      H5LTset_attribute_int( file_id, "eigenfunc", "iz", &izToDump, 1);
      H5LTset_attribute_int( file_id, "fdsolution", "iz", &izToDump, 1);
      H5LTset_attribute_int( file_id, "product", "iz", &izToDump, 1);
      H5LTset_attribute_int( file_id, "product", "mode", &modeToDump, 1);
      double returnIntegral = 0.5*integral*dx;
      H5LTset_attribute_double( file_id, "product", "integral", &returnIntegral, 1);
      H5Fclose(file_id);
      clog << "Debug: Integrand written to " << fname << endl;
    }
  #endif

  return integral*0.5*dx;
}
double CurvedWaveGuideFD::smoothedWG( double x, double z ) const
{
  static bool indicationWritten = false;
  if ( !indicationWritten )
  {
    clog << "\nUsing smoothed WG\n";
    indicationWritten = true;
  }

  double w = 0.1*width;
  double d = 2.0*R*x + z*z;
  double f1 = 1.0/( 1.0 + exp(-d/w) );

  double f2 = 1.0/( 1.0 + exp( (2*R*width-d)/w) );
  return 1.0 - (f1 - f2);
}

void CurvedWaveGuideFD::getXrayMatProp( double x, double z, double &delta, double &beta ) const
{
  delta = cladding->getDelta()*smoothedWG(x,z);
  beta = cladding->getBeta()*smoothedWG(x,z);
}
