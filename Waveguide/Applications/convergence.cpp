#include <iostream>
#include "cladding.hpp"
#include "curvedWaveGuide2D.hpp"
#include "crankNicholson.hpp"
#include "controlFile.hpp"
#include "straightWG2D.hpp"
#include <vector>
#include <H5Cpp.h>
#include <hdf5_hl.h>
#include <string>

using namespace std;

vector<double> errorRatio( const vector<double> &solutionAtTest )
{
  vector<double> error;
  for ( unsigned int i=1;i<solutionAtTest.size()-1;i++ )
  {
    error.push_back( abs(( solutionAtTest[i-1]-solutionAtTest[i] )/( solutionAtTest[i] - solutionAtTest[i+1] )) );
  }
  return error;
}

int main( int argc, char** argv )
{
  Cladding cladding1;
  double delta = 4.49E-5;
  double beta = 3.45E-6;
  cladding1.setRefractiveIndex(delta, beta);
  double width = 100.0; // Width of the waveguide in nm

  double zmin = 0.0;
  double zmax = 500.0E3; // 500 um
  double xmin = -0.2E3;
  double xmax = width + 0.2E3;

  const unsigned int Nx = 6;
  const unsigned int Nz = 6;
  double stepX[Nx] = {4.0, 2.0, 1.0, 0.5, 0.25, 0.125};
  double stepZ[Nz] = {400.0, 200.0, 100.0, 50.0, 25.0, 12.5};

  StraightWG2D wg1;
  double testZ = zmax/2.0;
  double testX = width/2.0;
  CrankNicholson solver;
  wg1.setWidth( width );
  wg1.setWaveLength( 0.1569 );
  wg1.setCladding( cladding1 );

  // Convergence when varying dx
  vector<double> solutionAtTestPointX;
  clog << "Varying stepsize in x-direction...\n";
  for ( unsigned int i=0;i<Nx;i++ )
  {
    clog << i+1 << " of " << Nx << endl;
    wg1.setTransverseDiscretization(xmin,xmax,stepX[i]);
    wg1.setLongitudinalDiscretization(zmin,zmax,stepZ[2]);
    wg1.setSolver(solver); // Need to set the solver to ensure after changing the discretization
    wg1.setBoundaryConditions();
    wg1.solve();
    solutionAtTestPointX.push_back( wg1.getIntensity(testX, testZ) );
  }

  vector<double> solutionAtTestPointZ;
  clog << "Varying stepsize in z-direction...\n";
  for ( unsigned int i=0;i<Nz;i++ )
  {
    clog << i+1 << " of " << Nz << endl;
    wg1.setTransverseDiscretization(xmin,xmax,stepX[2]);
    wg1.setLongitudinalDiscretization(zmin,zmax,stepZ[i]);
    wg1.setSolver(solver); // Need to set the solver to ensure after changing the discretization
    wg1.setBoundaryConditions();
    wg1.solve();
    solutionAtTestPointZ.push_back( wg1.getIntensity(testX, testZ) );
  }

  vector<double> errorX = errorRatio( solutionAtTestPointX );
  vector<double> errorZ = errorRatio( solutionAtTestPointZ );

  string fname("data/convergence.h5");
  hsize_t dim = errorX.size();
  int rank = 1;
  hid_t file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  H5LTmake_dataset( file_id, "errorRatioX", rank, &dim, H5T_NATIVE_DOUBLE, &errorX[0]);
  dim = Nx;
  H5LTmake_dataset( file_id, "stepsizeX", rank, &dim, H5T_NATIVE_DOUBLE, stepX);
  dim = errorZ.size();
  H5LTmake_dataset( file_id, "errorRatioZ", rank, &dim, H5T_NATIVE_DOUBLE, &errorZ[0]);
  dim = Nz;
  H5LTmake_dataset( file_id, "stepsizeZ", rank, &dim, H5T_NATIVE_DOUBLE, stepZ);
  H5Fclose(file_id);
  clog << "Convergence data writte to " << fname << endl;
  return 0;
}
