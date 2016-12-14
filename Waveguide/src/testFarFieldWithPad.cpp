#include "straightWG2D.hpp"
#include "planeWave.hpp"
#include "paraxialEquation.hpp"
#include "crankNicholson.hpp"
#include "cladding.hpp"
#include "controlFile.hpp"
#include <iostream>
#include <visa/visa.hpp>
#include <armadillo>
#include <chrono>
#include <thread>
#define KEEP_PLOT_FOR_SEC 10

using namespace std;
int main( int argc, char **argv )
{
  double wavelength = 0.1569;
  double width = 100.0;
  //wavelength = 0.5;
  double xmin = -width;
  double xmax = 2.0*width;
  double zmin = 0.0;
  double zmax = 100E3;
  unsigned int Nx = 1000;
  unsigned int Nz = 1000;
  double dx = (xmax-xmin)/Nx;
  double dz = (zmax-zmin)/Nz;

  StraightWG2D sim;
  Cladding cladding;
  try
  {
    sim.setWaveLength( wavelength );
    sim.setWidth( width );
    cladding.setRefractiveIndex(0.0,0.0);
    sim.setCladding(cladding);
    sim.setTransverseDiscretization( xmin, xmax, dx );
    sim.setLongitudinalDiscretization( zmin, zmax, dz );
    PlaneWave pw;
    pw.setWavelength( wavelength );
    pw.setAngleDeg( 0.05 );
    CrankNicholson solver;
    ParaxialEquation eq;
    solver.setEquation(eq);
    sim.setSolver( solver );
    sim.setBoundaryConditions( pw );
    clog << "Solving system...";
    sim.solve();
    clog << "done\n";
    sim.saveContour(false);
    ControlFile ctl("data/farFieldPadTest");
    sim.setFarFieldAngleRange( -0.1, 0.1 );
    sim.computeFarField( 1048576 );
    sim.save( ctl );
  }
  catch ( exception &exc )
  {
    cout << exc.what() << endl;
    return 1;
  }
  catch (...)
  {
    cout << "Unrecognized exception!\n";
    return 1;
  }

  return 0;
}
