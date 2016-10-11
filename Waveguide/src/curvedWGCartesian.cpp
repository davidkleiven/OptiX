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
  cdouble n(1.0-3.45E-5, 3.45E-6);
  cladding.setRefractiveIndex(n);
  double Rcurv = 40E6;
  double width = 100.0;
  double xmax = Rcurv + 4.0*width;
  double xmin = xmax - 5E3;
  double stepX = 1.0;
  double zmin = 0.0;
  double zmax = 500E3;
  double stepZ = 100.0;

  try
  {
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
    wg.solve();
    wg.save(fname);
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
  return 0;
}
