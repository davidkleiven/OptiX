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
//#define FAKE_RESULT_WITH_PURE_PHASE_SHIFT

using namespace std;

void generatePhaseObjectSolution( const map<string,double> &params, const arma::cx_mat &reference, arma::cx_mat &newSolution );

enum class Geometry_t{SPHERE, COATED_SPHERE};
int main( int argc, char **argv )
{
  cout << "--------------------\n";
  cout << "Modes:\n";
  cout << "0: Sphere\n";
  cout << "1: Coated sphere\n";
  cout << "--------------------\n";

  map<string,double> params;
  params["xmin_r"] = -4.0;            // Minimum x-position in units of the radius of the sphere
  params["xmax_r"] = 4.0;             // Maximum x-position in units of the radius of the sphere
  params["ymin_r"] = -4.0;            // Minimum y-position in units of the radius of the sphere
  params["ymax_r"] = 4.0;             // Maximum y-position in units of the radius of the sphere
  params["zmin_r"] = -2.0;            // Minimum z-position in units of the radius of the sphere
  params["zmax_r"] = 2.0;             // Maximum z-position in units of the radius of the sphere
  params["wavelength"] = 0.1569;      // Wavelength in nano meters
  params["radius"] = 2000.0;          // Radius of the sphere in nano meters
  params["Nt"] = 1024;                // Number of discretization points in the x and y direction (should 2^N)
  params["Nz"] = 512;                 // Number of discretization points in the z-direction (propagation direction)
  params["downSampleT"] = 8;          // Downsamplings factor in the x and y direction (1 --> no downsampling)
  params["downSampleZ"] = 4;          // Downsampling factor in the z-direction (1 --> no downsampling )
  params["savePlots"] = 0;            // Save the resulting plots to png file ( 0 --> no, otherwise yes )
  params["q_max_inverse_nm"] = 0.5;   // Maximum scattering wavevector to store for the farfield ( 1/nm )
  params["padLength"] = 2048;         // Length of the signal to use when computing the far field (should be 2^N)
  params["disableAbsorption"] = 0.0;  // Run without absorption (0 --> with absorption, otherwise run without absorption)
  params["mode"] = 0;
  params["coatThicknessIn_nm"] = 100.0; // Thickness of the coating in nano meters
  params["farFieldSaveSize"] = 512;     // Dimenstion of the farfield matrices to be saved
  params["intensityMinVis"] = 0.0;      // Minimum intensity for visualization
  params["intensityMaxVis"] = 1.1;      // Maximum intensity for visualization
  params["phaseMinVis"] = -1.58;        // Minimum phase for visualization
  params["phaseMaxVis"] = 1.58;         // Maximum phase for visualization

  pei::DialogBox dialog( params );
  dialog.show();

  Geometry_t geom;
  if ( static_cast<int>(params.at("mode")+0.5) == 0 )
  {
    geom = Geometry_t::SPHERE;
  }
  else if ( static_cast<int>(params.at("mode")+0.5) == 1 )
  {
    geom = Geometry_t::COATED_SPHERE;
  }

  double anglemax = params.at("q_max_inverse_nm")*params.at("wavelength")/(2.0*3.14159);
  anglemax *= (180.0/3.14159);

  // Define lengths in nm
  double r = params.at("radius");
  double xmin = params.at("xmin_r")*r;
  double xmax = params.at("xmax_r")*r;
  double zmin = params.at("zmin_r")*r;
  double zmax = params.at("zmax_r")*r;
  bool disableAbsorption = ( static_cast<int>(params.at("disableAbsorption")+0.5) != 0 );

  double dx = (xmax-xmin)/params.at("Nt");
  double dz = (zmax-zmin)/params.at("Nz");
  Point3D center;
  center.x = ( xmax+xmin)/2.0;
  center.y = center.x;
  center.z = (zmax+zmin)/2.0;

  post::Phase phase; // Only used for visualization
  post::ExitField ef;
  post::ExitIntensity ei;
  post::ExitPhase ep;
  post::FarField ff;

  unsigned int exportNx = 512;
  unsigned int exportNy = 512;
  ef.setExportDimensions( exportNx, exportNy );
  ei.setExportDimensions( exportNx, exportNy );
  ep.setExportDimensions( exportNx, exportNy );
  ff.setExportDimensions( params.at("farFieldSaveSize"), params.at("farFieldSaveSize") );

  Sphere *scatterer = NULL;
  Sphere sphere( center, r );
  CoatedSphere coated( center, r );
  coated.setThickness( params.at("coatThicknessIn_nm") );

  switch( geom )
  {
    case Geometry_t::SPHERE:
      scatterer = &sphere;
      break;
    case Geometry_t::COATED_SPHERE:
      scatterer = &coated;
      break;
  }

  if ( disableAbsorption )
  {
    scatterer->noAbsorption();
  }
  try
  {

    scatterer->setWaveLength( params.at("wavelength") );

    // Reference run
    scatterer->setMaterial( "Vacuum" );
    coated.setCoatingMaterial( "Vacuum" );

    scatterer->setTransverseDiscretization( xmin, xmax, dx, params.at("downSampleT") );
    scatterer->setVerticalDiscretization( xmin, xmax, dx );
    scatterer->setLongitudinalDiscretization( zmin, zmax, (zmax-zmin)/3.0, 1 );

    #ifdef FAKE_RESULT_WITH_PURE_PHASE_SHIFT
      FFT3DSolverDebug solver;
    #else
      FFTSolver3D solver;
    #endif

    PlaneWave pw;
    pw.setWavelength( params.at("wavelength") );
    pw.setDim( ParaxialSource::Dim_t::THREE_D );

    GaussianBeam gbeam;
    gbeam.setWavelength( params.at("wavelength") );
    gbeam.setDim( ParaxialSource::Dim_t::THREE_D );
    gbeam.setWaist( 400.0*r );

    scatterer->setSolver( solver );
    scatterer->setBoundaryConditions( gbeam );
    scatterer->solve();
    arma::cx_mat ref = solver.getLastSolution3D();
    clog << "Reference solution computed\n";

    // Compute real solution
    scatterer->reset();
    scatterer->setLongitudinalDiscretization( zmin, zmax, dz, params.at("downSampleZ") );
    scatterer->setMaterial( "SiO2" );
    //scatterer->setMaterial( "Au" );
    coated.setCoatingMaterial( "Au" );

    solver.visualizeRealSpace();
    solver.setIntensityMinMax( params.at("intensityMinVis"), params.at("intensityMaxVis") );
    solver.setPhaseMinMax( params.at("phaseMinVis"), params.at("phaseMaxVis") );
    if ( argc > 1 )
    {
      string imgname(argv[1]);
      solver.storeImages( imgname.c_str() );
    }

    scatterer->setSolver( solver );
    scatterer->setBoundaryConditions( gbeam );

    #ifndef FAKE_RESULT_WITH_PURE_PHASE_SHIFT
      clog << "Solving system...";
      scatterer->solve();
    clog << "done\n";
    #endif

    ControlFile ctl("data/sphere");

    #ifdef FAKE_RESULT_WITH_PURE_PHASE_SHIFT
      clog << "Using fake solution to test the far field routine!\n";
      arma::cx_mat fakeSolution;
      generatePhaseObjectSolution( params, ref, fakeSolution );
      solver.getLastSolution3D() = fakeSolution;
    #endif

    ff.setAngleRange( -anglemax, anglemax );
    ff.setPadLength( params.at("padLength")+0.5 );
    ff.setPadding( post::FarField::Pad_t::ZERO );
    ff.setReference( ref );
    *scatterer << ef << ei << ep << ff;
    scatterer->save( ctl );
    ctl.save();

    // Perform some plots of different projections
    visa::WindowHandler plots;
    plots.setLayout(2,3);
    typedef visa::Colormaps::Colormap_t cmap_t;
    plots.addPlot("IntensityXY");
    plots.addPlot("PhaseXY");
    plots.get("IntensityXY").setCmap( cmap_t::GREYSCALE );
    plots.get("PhaseXY").setCmap( cmap_t::GREYSCALE );
    unsigned int sliceNumber = params.at("Nz")/( 2.0*params.at("downSampleZ") );
    arma::mat solution = arma::abs( scatterer->getSolver().getSolution3D().slice(sliceNumber) );
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
    unsigned int row = params.at("Nt")/( 2.0*params.at("downSampleT") );
    solution = arma::abs( scatterer->getSolver().getSolution3D().tube( row, 0, row, params.at("Nt")/params.at("downSampleT")-1 ) );
    plots.get("IntensityYZ").fillVertexArray( solution );
    solution = arma::arg( scatterer->getSolver().getSolution3D().tube( row, 0, row, params.at("Nt")/params.at("downSampleT")-1 ) );
    plots.get("PhaseYZ").fillVertexArray( solution );

    // XZ plane
    plots.addPlot( "IntensityXZ" );
    plots.addPlot( "PhaseXZ" );
    plots.get("IntensityXZ").setCmap( cmap_t::GREYSCALE );
    plots.get("PhaseXZ").setCmap( cmap_t::GREYSCALE );
    unsigned int col = params.at("Nt")/(2.0*params.at("downSampleT"));
    solution = arma::abs( scatterer->getSolver().getSolution3D().tube( 0, col, params.at("Nt")/params.at("downSampleT")-1, col ) );
    plots.get("IntensityXZ").fillVertexArray( solution );
    solution = arma::arg( scatterer->getSolver().getSolution3D().tube( 0, col, params.at("Nt")/params.at("downSampleT")-1, col ) );
    plots.get("PhaseXZ").fillVertexArray( solution );
    plots.show();

    // Store all pictures
    if ( static_cast<int>(params.at("savePlots")) )
    {
      stringstream fname;
      fname << "Figures/sphere" << ctl.getUID() << ".jpg";
      plots.saveImg( fname.str().c_str() );
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

// Helper functions
void generatePhaseObjectSolution( const map<string,double> &params, const arma::cx_mat &reference, arma::cx_mat &newSolution )
{
  assert( reference.is_square() );
  double xmin = params.at("xmin_r")*params.at("radius");
  double xmax = params.at("xmax_r")*params.at("radius");
  double R = params.at("radius");
  double dx = (xmax-xmin)/reference.n_rows;
  double delta = 8.90652E-6;
  double k = 2.0*3.14159/params.at("wavelength");
  cdouble im( 0.0, 1.0 );
  newSolution.set_size( reference.n_rows, reference.n_cols );
  for ( unsigned int i=0;i<reference.n_cols;i++ )
  {
    double x = xmin + i*dx;
    for ( unsigned int j=0;j<reference.n_rows;j++ )
    {
      double y = xmin+j*dx;
      double r = sqrt( x*x+y*y );
      if ( r > R )
      {
        newSolution(j,i) = reference(j,i);
      }
      else
      {
        double phase = -2.0*delta*k*sqrt( R*R - r*r );
        phase = -2.0*delta*k*R;
        newSolution(j,i) = abs(reference(j,i))*exp(im*phase);
      }
      //newSolution(j,i) = 2.0*reference(j,i);
    }
  }
}
