#include <iostream>
#include "cladding.hpp"
#include "crankNicholson.hpp"
#include "controlFile.hpp"
#include "straightWG2D.hpp"
#include "coupledCurvedWG.hpp"
#include "curvedWGCylCrd.hpp"
#include "paraxialEquation.hpp"
#include "beamSplitter.hpp"
#include "planeWave.hpp"
#include <complex>
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>

using namespace std;

int main( int argc, char **argv )
{
  const double PI = acos(-1.0);
  Cladding cladding;
  double delta = 4.49E-5;
  double beta = 3.45E-6;

  // Values for SiN3
  delta = 1E-5;
  beta = 1E-7;
  cladding.setRefractiveIndex(delta, beta);
  double width = 100.0; // Width of the waveguide in nm

  double zmin = 0.0;
  double zmax = 800E3;
  double dz = 100.0;
  double splitStart = 100E3;
  double wgEnd = 800E3;
  double angleDeg = 0.015;
  double xmax = wgEnd*angleDeg*PI/180.0 + 1.5*width;
  double dx = 1.0;

  ControlFile ctl("data/beamSplitter");
  double wavelength = 0.1569;
  try
  {
    BeamSplitter wg;
    wg.setWaveLength( wavelength );
    wg.setCladding( cladding );
    wg.setTransverseDiscretization(-xmax,xmax,dx);
    wg.setLongitudinalDiscretization(zmin,zmax,dz);
    wg.setAngleWithZAxisDeg( angleDeg );
    wg.setSplitStart( splitStart );

    CrankNicholson solver;
    ParaxialEquation eq;
    PlaneWave pw;
    pw.setWavelength( wavelength );
    solver.setEquation( eq );
    wg.setSolver( solver );
    wg.setBoundaryConditions( pw );

    clog << "Solving equation... ";
    wg.solve();
    clog << "done\n";
    wg.extractWGBorders();
    wg.computeFarField( 524288 );
    clog << "Exporting results...\n";
    wg.save( ctl );
    ctl.save();
    clog << "Finished exporting\n";
  }
  catch ( exception &exc )
  {
    cerr << exc.what() << endl;
    return 1;
  }

  clog << "The simulation ended successfully\n";
  return 0;
}
