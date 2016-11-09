#include <iostream>
#include "cladding.hpp"
#include "curvedWaveGuide2D.hpp"
#include "crankNicholson.hpp"
#include "controlFile.hpp"
#include "straightWG2D.hpp"
#include "paraxialEquation.hpp"
#include "planeWave.hpp"
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

double compare( const WaveGuideFDSimulation &poorDisc, const WaveGuideFDSimulation &goodDisc )
{
  // Loop over the matrix with fewest entries
  unsigned int Nr = poorDisc.nodeNumberTransverse();
  unsigned int Nc = poorDisc.nodeNumberLongitudinal();
  double err1Sq = 0.0;
  double err2Sq = 0.0;
  double crossTerm = 0.0;
  for ( unsigned int ic=0;ic<Nc;ic++ )
  {
    for ( unsigned int ir=0;ir<Nr;ir++ )
    {
      double err1 = abs( poorDisc.getSolver().getSolution()(ir,ic) );
      double x = poorDisc.getX(ir);
      double z = poorDisc.getZ(ic);
      double err2 = sqrt( goodDisc.getIntensity(x,z) );
      err1Sq += err1*err1;
      err1Sq += err2*err2;
      crossTerm += 2.0*err1*err2;
    }
  }
  return sqrt( err1Sq+err2Sq - crossTerm);
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

  const unsigned int nElem = 7;
  int N[nElem] = {100, 200, 400, 800, 1200, 2400, 4800}; // Not unsigned, due to HDF5

  StraightWG2D wg1;
  StraightWG2D wg2;
  double testZ = zmax/2.0;
  double testX = width/2.0;
  ParaxialEquation eq1;
  ParaxialEquation eq2;
  CrankNicholson solver;
  CrankNicholson solver2;
  PlaneWave pw1;
  PlaneWave pw2;
  pw1.setWavelength( 0.1569 );
  pw2.setWavelength( 0.1569 );
  solver.setEquation(eq1);
  solver2.setEquation(eq2);
  wg1.setWidth( width );
  wg1.setWaveLength( 0.1569 );
  wg1.setCladding( cladding1 );
  wg2.setWidth( width );
  wg2.setWaveLength( 0.1569 );
  wg2.setCladding( cladding1 );

  // Convergence when varying dx
  vector<double> error;
  clog << "Varying stepsize in x-direction...\n";
  for ( unsigned int i=0;i<nElem-1;i++ )
  {
    clog << i+1 << " of " << nElem-1 << endl;
    double stepX1 = (xmax-xmin)/static_cast<double>(N[i]);
    double stepZ1 = (zmax-zmin)/static_cast<double>(N[i]);
    double stepX2 = (xmax-xmin)/static_cast<double>(N[i+1]);
    double stepZ2 = (zmax-zmin)/static_cast<double>(N[i+1]);
    wg1.setTransverseDiscretization(xmin,xmax,stepX1);
    wg1.setLongitudinalDiscretization(zmin,zmax,stepZ1);
    wg1.setSolver(solver); // Need to set the solver to ensure after changing the discretization
    wg1.setBoundaryConditions(pw1);
    wg1.solve();

    wg2.setTransverseDiscretization(xmin, xmax, stepX2);
    wg2.setLongitudinalDiscretization(zmin, zmax, stepZ2);
    wg2.setSolver(solver2);
    wg2.setBoundaryConditions(pw2);
    wg2.solve();
    error.push_back( compare(wg1,wg2) );
  }

  string fname("data/convergence.h5");
  hsize_t dim = error.size();
  int rank = 1;
  hid_t file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  H5LTmake_dataset( file_id, "error", rank, &dim, H5T_NATIVE_DOUBLE, &error[0]);
  H5LTset_attribute_double( file_id, "error", "xmin", &xmin, 1);
  H5LTset_attribute_double( file_id, "error", "xmax", &xmax, 1);
  H5LTset_attribute_double( file_id, "error", "zmin", &zmin, 1);
  H5LTset_attribute_double( file_id, "error", "zmax", &zmax, 1);
  dim = nElem;
  H5LTmake_dataset( file_id, "nodes", rank, &dim, H5T_NATIVE_INT, N);
  H5Fclose(file_id);
  clog << "Convergence data writte to " << fname << endl;
  return 0;
}
