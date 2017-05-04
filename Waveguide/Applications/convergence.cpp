#include <iostream>
#include "cladding.hpp"
#include "curvedWaveGuide2D.hpp"
#include "crankNicholson.hpp"
#include "straightWG2D.hpp"
#include "paraxialEquation.hpp"
#include "planeWave.hpp"
#include <visa/visa.hpp>
#include <vector>
#include <H5Cpp.h>
#include <hdf5_hl.h>
#include <string>
#include <chrono>
#include <thread>

using namespace std;
enum class ErrorNorm_t{L1, L2};
vector<double> errorRatio( const vector<double> &solutionAtTest )
{
  vector<double> error;
  for ( unsigned int i=1;i<solutionAtTest.size()-1;i++ )
  {
    error.push_back( abs(( solutionAtTest[i-1]-solutionAtTest[i] )/( solutionAtTest[i] - solutionAtTest[i+1] )) );
  }
  return error;
}

double compare( const WaveGuideFDSimulation &poorDisc, const WaveGuideFDSimulation &goodDisc, ErrorNorm_t norm )
{
  // Loop over the matrix with fewest entries
  unsigned int Nr = poorDisc.nodeNumberTransverse();
  unsigned int Nc = poorDisc.nodeNumberLongitudinal();
  double error = 0.0;
  double err2Sq = 0.0;
  double crossTerm = 0.0;

  for ( unsigned int ic=0;ic<Nc-1;ic++ )
  {
    for ( unsigned int ir=0;ir<Nr-1;ir++ )
    {
      double err1 = abs( poorDisc.getSolver().getSolution()(ir,ic) );
      double x = poorDisc.getX(ir);
      double z = poorDisc.getZ(ic);
      unsigned int ix, iz;
      goodDisc.closestIndex( x, z, ix, iz);
      double err2 = abs( goodDisc.getSolver().getSolution()(ix,iz) );

      switch ( norm )
      {
        case ErrorNorm_t::L1:
          error += abs(err1-err2);
          break;
        case ErrorNorm_t::L2:
          error += pow(err1-err2, 2);
          break;
      }
    }
  }

  switch ( norm )
  {
    case ErrorNorm_t::L1:
      return error/static_cast<double>(Nr*Nc);
    case ErrorNorm_t::L2:
      return sqrt( error/static_cast<double>(Nr*Nc) );
  }
}

int main( int argc, char** argv )
{
  bool visualizeResults = true;
  ErrorNorm_t norm = ErrorNorm_t::L1;
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
  visa::WindowHandler plots;

  if ( visualizeResults )
  {
    plots.addPlot("Active");
    plots.addPlot("Passive");
  }

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
    if (( i == 0 ) || ( i%2 == 1))
    {
      clog << "Wg1 uses: n = " << N[i+i%2] << endl;
      double stepX1 = (xmax-xmin)/static_cast<double>(N[i+i%2]);
      double stepZ1 = (zmax-zmin)/static_cast<double>(N[i+i%2]);
      wg1.setTransverseDiscretization(xmin,xmax,stepX1);
      wg1.setLongitudinalDiscretization(zmin,zmax,stepZ1);
      wg1.setSolver(solver); // Need to set the solver to ensure after changing the discretization
      wg1.setBoundaryConditions(pw1);
      wg1.solve();

      if ( visualizeResults )
      {
        arma::mat solutionCopy = arma::abs( wg1.getSolver().getSolution() );
        plots.get("Active").fillVertexArray( solutionCopy );
      }
    }

    if (( i == 0 ) || (i%2 == 0))
    {
      clog << "Wg2 uses: n = " << N[i+1] << endl;
      double stepX2 = (xmax-xmin)/static_cast<double>(N[i+1]);
      double stepZ2 = (zmax-zmin)/static_cast<double>(N[i+1]);

      wg2.setTransverseDiscretization(xmin, xmax, stepX2);
      wg2.setLongitudinalDiscretization(zmin, zmax, stepZ2);
      wg2.setSolver(solver2);
      wg2.setBoundaryConditions(pw2);
      wg2.solve();

      if ( visualizeResults )
      {
        arma::mat solutionCopy = arma::abs( wg2.getSolver().getSolution() );
        plots.get("Passive").fillVertexArray( solutionCopy );
      }
    }

    if (( i== 0 ) || (i%2 == 0))
    {
      clog << "Comparing using wg2 finer than wg1\n";
      error.push_back( compare(wg1, wg2, norm) );
    }
    else if ( i%2 == 1 )
    {
      clog << "Comparing using that wg1 finer than wg2\n";
      error.push_back( compare(wg2, wg1, norm) );
    }

    if ( visualizeResults )
    {
      // Very slow over SSH
      //for ( unsigned int i=0;i<10;i++ )
      {
        plots.show();
        //this_thread::sleep_for(chrono::milliseconds(300));
      }
    }
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
