#include "cylinder2D.hpp"
#include "planeWave.hpp"
#include "paraxialEquation.hpp"
#include "crankNicholson.hpp"
#include "controlFile.hpp"
#include "postProcessMod.hpp"
#include <iostream>
#include <visa/visa.hpp>
#include <armadillo>
#include <chrono>
#include <thread>
#define KEEP_PLOT_FOR_SEC 10

using namespace std;
int main( int argc, char **argv )
{
  // Define lengths in nm
  double x0 = 0.0;
  double radius = 2000.0;
  double z0 = 2.0*radius;
  double wavelength = 0.1569;
  //wavelength = 0.5;
  double xmin = -1.5*radius;
  double xmax = 1.5*radius;
  double zmin = 0.0;
  double zmax = 4.0*radius;
  unsigned int Nx = 1000;
  unsigned int Nz = 1000;
  double dx = (xmax-xmin)/Nx;
  double dz = (zmax-zmin)/Nz;

  post::Phase phase;
  Cylinder2D cylinderSim(x0,z0,radius);
  try
  {
    cylinderSim.setWaveLength( wavelength );
    cylinderSim.setMaterial( "SiO2" );
    cylinderSim.setTransverseDiscretization( xmin, xmax, dx );
    cylinderSim.setLongitudinalDiscretization( zmin, zmax, dz );
    PlaneWave pw;
    pw.setWavelength( wavelength );
    CrankNicholson solver;
    ParaxialEquation eq;
    solver.setEquation(eq);
    cylinderSim.setSolver( solver );
    cylinderSim.setBoundaryConditions( pw );
    clog << "Solving system...";
    cylinderSim.solve();
    clog << "done\n";
    ControlFile ctl("data/cylinder");

    post::FarField ff;
    ff.setAngleRange( -0.1, 0.1 );
    ff.setPadLength( 262144 );
    cylinderSim.save( ctl );

    visa::WindowHandler plots;
    plots.addPlot("Intensity");
    plots.addPlot("Phase");
    plots.get("Intensity").setColorMax(1.0);
    plots.get("Intensity").setColorMin(0.99);
    arma::mat solution = arma::abs( cylinderSim.getSolver().getSolution() );
    plots.get("Intensity").fillVertexArray( solution );
    plots.get("Phase").setColorMax(3.14159);
    plots.get("Phase").setColorMin(-3.14159);
    phase.result( solver, solution );
    plots.get("Phase").fillVertexArray(solution);

    plots.show();
    for ( unsigned int i=0;i<KEEP_PLOT_FOR_SEC;i++ )
    {
      plots.show();
      this_thread::sleep_for( chrono::seconds(1) );
      clog << "Closes in " << KEEP_PLOT_FOR_SEC-i << " seconds \r";
    }
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
