#include <iostream>
#include "cladding.hpp"
#include "curvedWaveGuide2D.hpp"
#include "crankNicholson.hpp"
#include "fftSolver2D.hpp"
#include "controlFile.hpp"
#include "straightWG2D.hpp"
#include "paraxialSource.hpp"
#include "planeWave.hpp"
#include "paraxialEquation.hpp"
#include "gaussianBeam.hpp"
#include "cylindricalParaxialEquation.hpp"
#include "curvedWGCylCrd.hpp"
#include "gaussianWG.hpp"
#include "linearRampWG.hpp"
#include "postProcessMod.hpp"
#include "curvedWGConfMap.hpp"
#include <map>
#include <complex>
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>
#include <visa/visa.hpp>
#include <chrono>
#include <thread>
#include <pei/dialogBox.hpp>
#define KEEP_PLOT_FOR_SEC 6
#define VISUALIZE_PATTERN

using namespace std;

typedef complex<double> cdouble;

enum class Source_t {PLANE, GAUSSIAN};
enum class WGProfile_t {STEP, GAUSSIAN, LINEAR_RAMP};
enum class Mode_t {STRAIGHT, CURVED, CURVED_CYLCRD, CURVED_CONF_MAP, FFT_CARTESIAN, NONE};

/** Function for the common setup for the wave guide */
void commonSetup( CurvedWaveGuideFD &wg, const map<string, double> &params )
{
  wg.setRadiusOfCurvature( params.at("radius")*1E6 );
  wg.setWaveguideLength( params.at("wglength") );
  wg.setWidth( params.at("width") );
  wg.setWaveLength( params.at("wavelength") );
  double stepX = (params.at("xmax") - params.at("xmin"))/params.at("Nx");
  double stepZ = (params.at("zmax") - params.at("zmin"))/params.at("Nz");
  wg.setTransverseDiscretization( params.at("xmin"), params.at("xmax"), stepX, params.at("downSamplingX") );
  wg.setLongitudinalDiscretization( params.at("zmin"), params.at("zmax"), stepZ, params.at("downSamplingZ") );
}

/** Setup the plane wave function */
void setupPW( PlaneWave &pw, const map<string,double> &params )
{
  pw.setWavelength( params.at("wavelength") );
  pw.setAngleDeg( params.at("incAngle") );
}

/** Main function */
typedef visa::Colormaps::Colormap_t cmap_t;

int main( int argc, char **argv )
{
  map<string,double> params;
  params["radius"] = 40.0;             // In mm
  params["wavelength"] = 0.1569;       // In nm
  params["width"] = 100.0;             // In nm
  params["Nx"] = 2000;
  params["Nz"] = 6000;
  params["zmin"] = 0.0;
  params["zmax"] = 400E3;              // in nm
  params["ffAngleMax"] = 1.0;          // In deg
  params["ffAngleMin"] = -1.0;         // In deg
  params["incAngle"] = 0.0;            // In deg
  params["xmin"] = 0.0;
  params["xmax"] = 0.0;
  params["downSamplingX"] = 10.0;
  params["downSamplingZ"] = 10.0;
  params["wglength"] = 1.1*params["zmax"];
  params["padLength"] = 131072;
  params["mode"] = 2.0;
  params["borderTracker"] = 0;
  Mode_t mode = Mode_t::NONE;

  cout << "Mode description:\n";
  cout << "    1: Straight\n";
  cout << "    2: Curved waveguide in cylindrical coordinates\n";
  cout << "    3: Curved waveguide in cartesian coordinates\n";
  cout << "    4: Curved waveguide using conformal map\n";
  cout << "    5: Curved waveguide in cartesian coordinates using FFT solver\n";
  cout << "brdtrack: Use the border tracker. Only have effect if using cartesian coordinates\n";
  cout << flush;

  pei::DialogBox box( params );
  box.show();

  bool useBorderTracker = static_cast<int>( params["borderTracker"]+0.5 );
  unsigned int intMode = static_cast<unsigned int>( params["mode"]+0.5 );
  if ( intMode == 1 ) mode = Mode_t::STRAIGHT;
  else if ( intMode == 2 ) mode = Mode_t::CURVED_CYLCRD;
  else if ( intMode == 3 ) mode = Mode_t::CURVED;
  else if ( intMode == 4 ) mode = Mode_t::CURVED_CONF_MAP;
  else if ( intMode == 5 ) mode = Mode_t::FFT_CARTESIAN;

  Cladding cladding;
  double delta = 4.14E-5; // Default value
  double beta = 3.45E-6;
  cladding.setRefractiveIndex(delta, beta);

  // Set post processing options
  post::Intensity amplitude;
  post::Phase phase;
  post::ExitField ef;
  post::ExitIntensity ei;
  post::ExitPhase ep;
  post::FarField ff;
  ff.setPadLength( params["padLength"] );
  ff.setAngleRange( params["ffAngleMin"], params["ffAngleMax"] );
  ControlFile ctl("data/singleCurvedWG"); // File for all parameters and settings
  PlaneWave pw;
  setupPW( pw, params );
  CrankNicholson solver;
  solver.setBoundaryCondition( Solver2D::BC_t::TRANSPARENT );
  arma::mat intensity; // Matrix used for visualization
  double xMargin = 0.01E3;

  try
  {
    switch ( mode )
    {
      case Mode_t::STRAIGHT:
      {
        clog << "Running straight waveguide\n";
        params["xmin"] = -params["width"];
        params["xmax"] = 2.0*params["width"];
        StraightWG2D wg;
        commonSetup( wg, params );
        ParaxialEquation eq;
        solver.setEquation( eq );
        wg.setCladding( cladding );
        wg.setSolver( solver );
        wg.setBoundaryConditions( pw );

        // Set post processing options
        wg << amplitude << phase << ef << ei << ep << ff;
        clog << "Solving linear system... ";
        wg.solve();
        clog << "done\n";

        clog << "Exporting results...\n";
        wg.save( ctl );
        ctl.save();
        clog << "Finished exporting\n";
        intensity = arma::abs( wg.getSolver().getSolution() );
        break;
      }
      case Mode_t::CURVED_CYLCRD:
      {
        clog << "Solving curved waveguide in cylindrical coordinates\n";
        params["xmin"] = -params["width"];
        params["xmax"] = 2.0*params["width"];
        params["zmax"] /= ( params["radius"]*1E6 ); // Convert to radians
        CurvedWGCylCrd wg;
        commonSetup( wg, params );
        CylindricalParaxialEquation eq;
        eq.setRadiusOfCurvature( params["radius"]*1E6 );
        solver.setEquation( eq );
        wg.setCladding( cladding );
        wg.setSolver( solver );
        wg.setBoundaryConditions( pw );
        wg << amplitude << phase << ef << ei << ep << ff;
        clog << "Solving equation...";
        wg.solve();
        clog << " done\n";
        clog << "Exporting results...\n";
        wg.save( ctl );
        ctl.save();
        clog << "Finished exporting\n";
        intensity = arma::abs( wg.getSolver().getSolution() );
        break;
      }
      case Mode_t::CURVED:
      {
        clog << "Curved waveguide in cartesian coordinates\n";
        CurvedWaveGuideFD wg;
        params["xmax"] = params["width"] + xMargin;
        params["xmin"] = -0.5*params["zmax"]*params["zmax"]/(params["radius"]*1E6) - xMargin;
        commonSetup( wg, params );
        wg.setCladding( cladding );
        ParaxialEquation eq;
        solver.setEquation( eq ),
        wg.setSolver( solver );
        wg.setBoundaryConditions( pw );
        wg << amplitude << phase << ef << ei << ep << ff;
        if ( useBorderTracker ) wg.useBorderTracker();

        clog << "Solving system...";
        wg.solve();
        clog << " done\n";
        clog << "Exporting results...\n";
        if (!useBorderTracker) wg.extractWGBorders();
        wg.save( ctl );
        ctl.save();
        clog << "Finished exporting results\n";
        intensity = arma::abs( wg.getSolver().getSolution() );
        break;
      }
      case Mode_t::CURVED_CONF_MAP:
      {
        clog << "Curved waveguide in uv-plane\n";
        CurvedWGConfMap wg;
        params["xmax"] = params["width"];
        params["xmin"] = -2.0*params["width"];
        commonSetup( wg, params );
        wg.setCladding( cladding );
        ParaxialEquation eq;
        solver.setEquation( eq );
        wg.setSolver( solver );
        wg.setBoundaryConditions( pw );
        wg << amplitude << phase << ef << ei << ep << ff;

        clog << "Solving system...";
        wg.solve();
        clog << " done\n";
        clog << "Exporting results...\n";
        if (!useBorderTracker) wg.extractWGBorders();
        wg.save( ctl );
        ctl.save();
        clog << "Finished exporting results\n";
        intensity = arma::abs( wg.getSolver().getSolution() );
        break;
      }
      case Mode_t::FFT_CARTESIAN:
      {
        clog << "Curved waveguide in cartesian coordinates\n";
        CurvedWaveGuideFD wg;
        params["xmax"] = params["width"] + xMargin;
        params["xmin"] = -0.5*params["zmax"]*params["zmax"]/(params["radius"]*1E6) - xMargin;
        commonSetup( wg, params );
        wg.setCladding( cladding );
        FFTSolver2D fftSolver;
        wg.setSolver( fftSolver );
        wg.setBoundaryConditions( pw );
        wg << amplitude << phase << ef << ei << ep << ff;

        clog << "Solving system...";
        wg.solve();
        clog << " done\n";
        clog << "Exporting results...\n";
        if (!useBorderTracker) wg.extractWGBorders();
        wg.save( ctl );
        ctl.save();
        clog << "Finished exporting results\n";
        intensity = arma::abs( wg.getSolver().getSolution() );
        break;
      }
        default:
          clog << "Mode not recognized\n";
    }

    #ifdef VISUALIZE_PATTERN
      clog << "Visualizing waveguide intensity\n";
      visa::WindowHandler plots;
      plots.addPlot("Intensity");
      plots.get("Intensity").setCmap( cmap_t::NIPY_SPECTRAL );
      plots.get("Intensity").fillVertexArray( intensity );
      plots.show();
      for ( unsigned int i=0;i<KEEP_PLOT_FOR_SEC;i++ )
      {
        plots.show();
        clog << "Closes in " << KEEP_PLOT_FOR_SEC-i <<  "seconds \r";
        this_thread::sleep_for( chrono::seconds(1) );
      }
    #endif
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
  return 0;
}
