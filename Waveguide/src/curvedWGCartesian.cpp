#include <iostream>
#include "cladding.hpp"
#include "curvedWaveGuide2D.hpp"
#include "crankNicholson.hpp"
#include "controlFile.hpp"
#include <complex>
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>

using namespace std;

typedef complex<double> cdouble;

int main( int argc, char **argv )
{

  Cladding cladding;
  ControlFile ctl("data/singleCurvedWG"); // File for all parameters and settings
  //double delta = 4.14E-5;
  double delta = 4.49E-5;
  double beta = 3.45E-6;
  cladding.setRefractiveIndex(delta, beta);
  double Rcurv = 40E6;
  double width = 100.0;
  //double xmax = Rcurv + 1E3;//4.0*width;
  //double xmin = xmax - 2E3;
  double xmax = width+0.5E3;//4.0*width;
  double xmin = xmax-5.0E3;
  double stepX = 1.0;
  double zmin = 0.0;
  double zmax = 500E3;
  double stepZ = 100.0;

  try
  {
    clog << "Initializing simulation...";

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
    wg.computeTransmission( 5E4 ); // Compute transmission every 50 um
    clog << "done\n";
    clog << "Exporting results...\n";
    wg.save( ctl);
    ctl.save();
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
