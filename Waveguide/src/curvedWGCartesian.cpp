#include <iostream>
#include "cladding.hpp"
#include "curvedWaveGuide2D.hpp"
#include "crankNicholson.hpp"
#include <complex>
#include <stdexcept>

using namespace std;

typedef complex<double> cdouble;

int main( int argc, char **argv )
{
  Cladding cladding;
  double delta = 3.45E-5;
  double beta = 3.45E-6;
  cladding.setRefractiveIndex(delta, beta);
  double Rcurv = 40E6;
  double width = 100.0;
  double xmax = Rcurv + 4.0*width;
  double xmin = xmax - 5E3;
  double stepX = 10.0;
  double zmin = 0.0;
  double zmax = 5E3;
  double stepZ = 1000.0;

  try
  {
    clog << "Initializing simulation...";
    string fname( "data/intensity2D_FD" );
    CurvedWaveGuideFD wg;
    wg.setRadiusOfCurvature( 40E6 );
    wg.setWidth( 100.0 );
    wg.setWaveLength( 0.1569 );
    wg.setCladding(cladding);
    wg.setTransverseDiscretization(xmin,xmax,stepX);
    wg.setLongitudinalDiscretization(zmin,zmax,stepZ);
    CrankNicholson solver;
    wg.setSolver(solver);
    wg.setBoundaryConditions();
    clog << " done\n";
    clog << "Solving linear system... ";
    wg.solve();
    clog << "done\n";
    clog << "Exporting results...\n";
    wg.save(fname);
    clog << "Finished exporting\n";

  }
  catch ( exception &exc )
  {
    cerr << exc.what() << endl;
    return 1;
  }
  catch (...)
  {
    cerr << "An unrecognized exception occured!\n";
  }

  clog << "The simulation ended successfully\n";
  return 0;
}
