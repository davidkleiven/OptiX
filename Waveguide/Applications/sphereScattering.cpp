#include "sphere.hpp"
#include "planeWave.hpp"
#include "paraxialEquation.hpp"
#include "crankNicholson.hpp"
#include "controlFile.hpp"
#include "postProcessMod.hpp"
#include "fftSolver3D.hpp"
#include <iostream>
#include <visa/visa.hpp>
#include <pei/dialogBox.hpp>
#include <armadillo>
#include <chrono>
#include <thread>
#include <map>
#define KEEP_PLOT_FOR_SEC 10

using namespace std;
int main( int argc, char **argv )
{
  map<string,double> params;
  params["xmin_r"] = -4.0;
  params["xmax_r"] = 4.0;
  params["ymin_r"] = -4.0;
  params["ymax_r"] = 4.0;
  params["zmin_r"] = -2.0;
  params["zmax_r"] = 2.0;
  params["wavelength"] = 0.1569;
  params["radius"] = 2000.0;
  params["Nt"] = 1024;
  params["Nz"] = 512;
  params["downSampleT"] = 8;
  params["downSampleZ"] = 4;

  pei::DialogBox dialog( params );
  dialog.show();

  // Define lengths in nm
  double r = params.at("radius");
  double xmin = params.at("xmin_r")*r;
  double xmax = params.at("xmax_r")*r;
  double zmin = params.at("zmin_r")*r;
  double zmax = params.at("zmax_r")*r;

  double dx = (xmax-xmin)/params.at("Nt");
  double dz = (zmax-zmin)/params.at("Nz");
  Point3D center;
  center.x = ( xmax+xmin)/2.0;
  center.y = center.x;
  center.z = (zmax+zmin)/2.0;

  post::Phase phase;
  post::Intensity amplitude;
  post::ExitField ef;
  post::ExitIntensity ei;
  post::ExitPhase ep;
  post::FarField ff;

  Sphere sphere( center, r );
  try
  {
    sphere.setWaveLength( params.at("wavelength") );
    sphere.setMaterial( "SiO2" );
    sphere.setTransverseDiscretization( xmin, xmax, dx, 1 );
    sphere.setLongitudinalDiscretization( zmin, zmax, dz, 1 );
    PlaneWave pw;
    pw.setWavelength( params.at("wavelength") );
    pw.setDim( ParaxialSource::Dim_t::THREE_D );

    FFTSolver3D solver;
    sphere.setSolver( solver );
    sphere.setBoundaryConditions( pw );
    clog << "Solving system...";
    sphere.solve();
    clog << "done\n";
    ControlFile ctl("data/sphere");

    ff.setAngleRange( -0.1, 0.1 );
    ff.setPadLength( 65536 );
    sphere << amplitude << phase << ef << ei << ep;
    sphere.save( ctl );
    ctl.save();

    visa::WindowHandler plots;
    plots.addPlot("Intensity");
    plots.addPlot("Phase");
    unsigned int sliceNumber = params.at("Nz")/2;
    arma::mat solution = arma::abs( sphere.getSolver().getSolution3D().slice(sliceNumber) );
    plots.get("Intensity").fillVertexArray( solution );
    arma::cube phaseField;
    phase.result( solver, phaseField );
    solution = phaseField.slice(sliceNumber);
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
