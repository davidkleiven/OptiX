#include <iostream>
#include "cladding.hpp"
#include "curvedWaveGuide2D.hpp"
#include "crankNicholson.hpp"
#include <complex>
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>

#define UID_DIGITS 6

using namespace std;
unsigned int UID_MAX = pow(10,UID_DIGITS);

typedef complex<double> cdouble;

int main( int argc, char **argv )
{
  srand(time(NULL));
  unsigned int uid = rand()%UID_MAX;

  Cladding cladding;
  double delta = 4.14E-5;
  double beta = 3.45E-6;
  //double delta = 0.0;
  //double beta = 0.0;
  cladding.setRefractiveIndex(delta, beta);
  double Rcurv = 40E6;
  double width = 100.0;
  //double xmax = Rcurv + 1E3;//4.0*width;
  //double xmin = xmax - 2E3;
  double xmax = width+1E3;//4.0*width;
  double xmin = -4E3;
  double stepX = 3.0;
  double zmin = 0.0;
  double zmax = 500E3;
  double stepZ = 3.0;

  try
  {
    clog << "UID for the current run " << uid << endl;
    clog << "Initializing simulation...";
    string fname( "data/intensity2D_FD" );
    stringstream ss;
    ss << fname;
    ss << uid;
    fname = ss.str();

    CurvedWaveGuideFD wg;
    wg.setRadiusOfCurvature( Rcurv );
    wg.setWidth( width );
    wg.setWaveLength( 0.1569 );
    wg.setCladding( cladding );
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