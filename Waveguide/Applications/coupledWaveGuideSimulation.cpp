#include <iostream>
#include "cladding.hpp"
#include "curvedWaveGuide2D.hpp"
#include "crankNicholson.hpp"
#include "controlFile.hpp"
#include "straightWG2D.hpp"
#include "coupledCurvedWG.hpp"
#include "curvedWGCylCrd.hpp"
#include "cylindricalParaxialEquation.hpp"
#include "planeWave.hpp"
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
  bool useCyl = true;
  double R = 80.0;

  // Parameters for running a sweep over radii of curvature
  double LzOverR = 0.01; // max(z)/R << 1 is a requirement
  double xMarginAboveAndBelow = 0.5E3; // In nanometers = 0.5 um
  unsigned int Nz = 5000; // Number of discretization points in x and z direction
  unsigned int nPointsTransmission = 200;

  Cladding cladding;
  double delta = 4.49E-5;
  double beta = 3.45E-6;
  cladding.setRefractiveIndex(delta, beta);
  double width = 100.0; // Width of the waveguide in nm

  // Start sweep

  double Rcurv = R*1E6;
  double zmin = 0.0;
  double zmax = Rcurv*LzOverR;
  double xmax = width+xMarginAboveAndBelow;
  double xmin = -0.5*zmax*LzOverR-xMarginAboveAndBelow;

  if ( useCyl )
  {
    zmin = 0.0;
    zmax = 0.01;
    xmin = -2.0*width;
    xmax = 6.0*width;
  }

  double stepX = (xmax-xmin)/static_cast<double>(Nz);
  double stepZ = (zmax-zmin)/static_cast<double>(Nz);
  stepX = stepX > 1.0 ? 1.0:stepX;
  stepZ = stepZ > 100.0 ? 100.0:stepZ;

  ControlFile ctl("data/coupledCurvedWG"); // File for all parameters and settings
  CoupledCurvedWG wg(CoupledCurvedWG::Coordinate_t::CYLINDRICAL);
  try
  {
    clog << "Initializing simulation...";

    wg.getWg1().setRadiusOfCurvature( Rcurv );
    wg.getWg2().setRadiusOfCurvature( Rcurv );
    wg.getWg1().setWidth( width );
    wg.getWg2().setWidth( width );

    double wavelength = 0.1569;
    wg.setWaveLength( wavelength );
    wg.setCladding( cladding );
    wg.setTransverseDiscretization(xmin,xmax,stepX);
    wg.setLongitudinalDiscretization(zmin,zmax,stepZ);
    wg.setSeparation( 1.0 );
    wg.setStartCoupler( 50E3/Rcurv );

    CrankNicholson solver;
    CylindricalParaxialEquation eq;
    PlaneWave pw;
    pw.setWavelength( wavelength );
    eq.setRadiusOfCurvature( Rcurv );
    solver.setEquation( eq );
    wg.setSolver(solver);
    wg.setBoundaryConditions( pw );
    clog << " done\n";
    clog << "Solving linear system... ";
    wg.solve();
    clog << "done\n";
    wg.extractWGBorders();
    /*
    clog << "Computing transmission... ";
    wg->computeTransmission( (zmax-zmin)/static_cast<double>(nPointsTransmission) );
    clog << "done\n";
    */
    clog << "Exporting results...\n";
    wg.save( ctl );
    //wg->saveTransmission( ctl );
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
    return 1;
  }

  clog << "The simulation ended successfully\n";
  return 0;
}
