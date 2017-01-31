#include "sphere.hpp"
#include "planeWave.hpp"
#include "paraxialEquation.hpp"
#include "gaussianBeam.hpp"
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
#define KEEP_PLOT_FOR_SEC 3

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
  params["savePlots"] = 0;
  params["q_max_inverse_nm"] = 0.5;
  params["padLength"] = 32768;

  pei::DialogBox dialog( params );
  dialog.show();

  double anglemax = params.at("q_max_inverse_nm")*params.at("wavelength")/(2.0*3.14159);
  anglemax *= (180.0/3.14159);

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

    // Reference run
    sphere.setMaterial( "SiO2" );
    sphere.setTransverseDiscretization( xmin, xmax, dx, 1 );
    sphere.setLongitudinalDiscretization( zmin, zmax, (zmax-zmin)/3.0, 1 );
    FFTSolver3D solver;

    PlaneWave pw;
    pw.setWavelength( params.at("wavelength") );
    pw.setDim( ParaxialSource::Dim_t::THREE_D );

    GaussianBeam gbeam;
    gbeam.setWavelength( params.at("wavelength") );
    gbeam.setDim( ParaxialSource::Dim_t::THREE_D );
    gbeam.setWaist( 400.0*r );

    sphere.setSolver( solver );
    sphere.setBoundaryConditions( gbeam );
    sphere.solve();
    arma::cx_mat ref = solver.getLastSolution3D();
    clog << "Reference solution computed\n";

    // Compute real solution
    sphere.reset();
    sphere.setLongitudinalDiscretization( zmin, zmax, dz, 1 );
    sphere.setMaterial( "SiO2" );

    solver.visualizeRealSpace();
    //solver.visualizeFourierSpace();
    solver.setIntensityMinMax( 0.0, 1.0 );
    solver.setPhaseMinMax( 0.0, 1.5 );

    sphere.setSolver( solver );
    sphere.setBoundaryConditions( gbeam );
    clog << "Solving system...";
    sphere.solve();
    clog << "done\n";
    ControlFile ctl("data/sphere");

    ff.setAngleRange( -anglemax, anglemax );
    ff.setPadLength( params.at("padLength")+0.5 );
    ff.setPadding( post::FarField::Pad_t::ZERO );
    ff.setReference( ref );
    sphere << ef << ei << ep << ff;
    sphere.save( ctl );
    ctl.save();

    // Perform some plots of different projections
    visa::WindowHandler plots;
    typedef visa::Colormaps::Colormap_t cmap_t;
    plots.addPlot("IntensityXY");
    plots.addPlot("PhaseXY");
    plots.get("IntensityXY").setCmap( cmap_t::GREYSCALE );
    plots.get("PhaseXY").setCmap( cmap_t::GREYSCALE );
    unsigned int sliceNumber = params.at("Nz")/2;
    arma::mat solution = arma::abs( sphere.getSolver().getSolution3D().slice(sliceNumber) );
    plots.get("IntensityXY").fillVertexArray( solution );
    arma::cube phaseField;
    phase.result( solver, phaseField );
    solution = phaseField.slice(sliceNumber);
    plots.get("PhaseXY").fillVertexArray(solution);

    // YZ plane
    plots.addPlot( "IntensityYZ" );
    plots.addPlot( "PhaseYZ" );
    plots.get("IntensityYZ").setCmap( cmap_t::GREYSCALE );
    plots.get("PhaseYZ").setCmap( cmap_t::GREYSCALE );
    unsigned int row = params.at("Nt")/2;
    solution = arma::abs( sphere.getSolver().getSolution3D().tube( row, 0, row, params.at("Nt") ) );
    plots.get("IntensityYZ").fillVertexArray( solution );
    solution = arma::arg( sphere.getSolver().getSolution3D().tube( row, 0, row, params.at("Nt") ) );
    plots.get("PhaseYZ").fillVertexArray( solution );

    // XZ plane
    plots.addPlot( "IntensityXZ" );
    plots.addPlot( "PhaseXZ" );
    plots.get("IntensityXZ").setCmap( cmap_t::GREYSCALE );
    plots.get("PhaseXZ").setCmap( cmap_t::GREYSCALE );
    unsigned int col = params.at("Nt")/2;
    solution = arma::abs( sphere.getSolver().getSolution3D().tube( 0, col, params.at("Nt"), col ) );
    plots.get("IntensityXZ").fillVertexArray( solution );
    solution = arma::arg( sphere.getSolver().getSolution3D().tube( 0, col, params.at("Nt"), col ) );
    plots.get("PhaseXZ").fillVertexArray( solution );

    // Store all pictures
    if ( static_cast<int>(params.at("savePlots")) )
    {
      for ( unsigned int i=0;i<plots.nPlots();i++ )
      {
        stringstream fname;
        fname << "Figures/sphere" << plots.get(i).getName() << ctl.getUID() << ".jpg";
        sf::Image img = plots.get(i).capture();
        img.saveToFile( fname.str() );
        clog << "Image " << fname.str() << " saved...\n";
      }
    }

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
