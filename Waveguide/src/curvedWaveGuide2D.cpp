#include "curvedWaveGuide2D.hpp"
#include <cassert>
#include <cmath>
#include <vector>
#include "cladding.hpp"
#include "solver2D.hpp"
#include <cmath>
#include "controlFile.hpp"
#include <H5Cpp.h>
#include <hdf5_hl.h>

using namespace std;

void CurvedWaveGuideFD::getXrayMatProp( double x, double z, double &delta, double &beta) const
{
  // Some assertions for debugging
  assert( cladding != NULL );
  assert ( x >= xDisc->min );
  assert ( x <= xDisc->max );
  assert ( z >= zDisc->min );
  assert ( z <= zDisc->max );
  if ( isInsideGuide( x, z) )
  {
    beta = 0.0;
    delta = 0.0;
    return;
  }
  delta = cladding->getDelta();
  beta = cladding->getBeta();
}

bool CurvedWaveGuideFD::isInsideGuide( double x, double z ) const
{
//  double d = sqrt(x*x + z*z);
  return (2.0*x*R+z*z > 0.0 ) && (2.0*x*R+z*z < 2.0*width*R);
  //return (pow(R+x,2)+z*z > pow(R,2)) && (pow(R+x,2)+z*z < pow(R+width,2));
  //return ( d < R+width) && (d > R );
}

void CurvedWaveGuideFD::setBoundaryConditions()
{
  unsigned int Nx = nodeNumberTransverse();
  vector<cdouble> values(Nx, 1.0);
  solver->setLeftBC(&values[0]);
}

void CurvedWaveGuideFD::fillInfo( Json::Value &obj ) const
{
  obj["RadiusOfCurvature"] = R;
  obj["Width"] = width;
}

cdouble CurvedWaveGuideFD::transverseBC( double z ) const
{
  double delta = cladding->getDelta();
  double beta = cladding->getBeta();
  cdouble im(0.0,1.0);
  return exp(-beta*wavenumber*z)*exp(-im*delta*wavenumber*z);
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
  double z = step;
  while ( z < zDisc->max )
  {
    double xWgStart = -0.5*z*z/R;
    closestIndex( xWgStart, z, wgStart, zIndx );
    closestIndex( xWgStart+width, z, wgEnd, zIndx );
    double intensity = trapezoidalIntegrateIntensityZ( zIndx, wgStart, wgEnd );
    transmission.push_back( intensity/intensityAtZero );
    z += step;
  }
  stepWhenComputingTransmission = step; // Save for later
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
