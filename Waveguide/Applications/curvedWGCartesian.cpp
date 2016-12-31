#include <iostream>
#include "cladding.hpp"
#include "curvedWaveGuide2D.hpp"
#include "crankNicholson.hpp"
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
#define KEEP_PLOT_FOR_SEC 6
#define VISUALIZE_PATTERN

using namespace std;

typedef complex<double> cdouble;

enum class Source_t {PLANE, GAUSSIAN};
enum class WGProfile_t {STEP, GAUSSIAN, LINEAR_RAMP};
enum class Mode_t {STRAIGHT, CURVED, CURVED_CYLCRD, NONE};

/** Function for the common setup for the wave guide */
void commonSetup( CurvedWaveGuideFD &wg, const map<string, double> &params )
{
  wg.setRadiusOfCurvature( params.at("radius")*1E6 );
  wg.setWaveguideLength( params.at("zmax") );
  wg.setWidth( params.at("width") );
  wg.setWaveLength( params.at("wavelength") );
  double stepX = (params.at("xmax") - params.at("xmin"))/params.at("Nx");
  double stepZ = (params.at("zmax") - params.at("zmin"))/params.at("Nz");
  wg.setTransverseDiscretization( params.at("xmin"), params.at("xmax"), stepX);
  wg.setLongitudinalDiscretization( params.at("zmin"), params.at("zmax"), stepZ);
}

/** Setup the plane wave function */
void setupPW( PlaneWave &pw, const map<string,double> &params )
{
  pw.setWavelength( params.at("wavelength") );
  pw.setAngleDeg( params.at("incAngle") );
}

/** Main function */
int main( int argc, char **argv )
{
  map<string,double> params;
  params["radius"] = 40.0;
  params["wavelength"] = 0.1569;
  params["width"] = 100.0;
  params["Nx"] = 2000;
  params["Nz"] = 6000;
  params["zmin"] = 0.0;
  params["zmax"] = 400.0;
  params["ffAngleMax"] = 1.0;
  params["ffAngleMin"] = -1.0;
  params["incAngle"] = 0.0;
  params["xmin"] = 0.0;
  params["xmax"] = 0.0;
  params["downSamplingX"] = 1.0;
  params["downSamplingZ"] = 1.0;
  Mode_t mode = Mode_t::NONE;

  bool useBorderTracker = false;

  /*********** PARSE COMMANDLINE ARGUMENTS ************************************/
  for ( unsigned int i=1;i<argc; i++ )
  {
    string arg(argv[i]);
    if ( arg.find("--help") != string::npos )
    {
      cout << "Usage: ./curvedWG.out [--help, --mode=<run number> --brdtrack]\n";
      cout << "help: Print this message\n";
      cout << "mode:\n";
      cout << "    1: Straight\n";
      cout << "    2: Curved waveguide in cylindrical coordinates\n";
      cout << "    3: Curved waveguide in cartesian coordinates\n";
      cout << "brdtrack: Use the border tracker. Only have effect if using cartesian coordinates\n";
      return 0;
    }
    else if ( arg.find("--mode=") != string::npos )
    {
      stringstream ss;
      unsigned int intMode;
      ss << arg.substr(7);
      ss >> intMode;
      if ( intMode == 1 ) mode = Mode_t::STRAIGHT;
      else if ( intMode == 2 ) mode = Mode_t::CURVED_CYLCRD;
      else if ( intMode == 3 ) mode = Mode_t::CURVED;
    }
    else if ( arg.find("--brdtrack") != string::npos )
    {
      useBorderTracker = true;
    }
    else
    {
      cout << "Unknown argument " << arg << endl;
      return 0;
    }
  }
  /****************** END COMMANDLINE ARGUMENTS *******************************/

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
  arma::mat intensity; // Matrix used for visualization
  double xMargin = 0.01E3;

  try
  {
    clog << "Initializing simulation...";
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
        CurvedWGCylCrd wg;
        commonSetup( wg, params );
        CylindricalParaxialEquation eq;
        eq.setRadiusOfCurvature( params["radius"]*1E6 );
        solver.setEquation( eq );

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
        params["xmax"] = 2.0*params["width"];
        params["xmin"] = -0.5*params["zmax"]*params["zmax"]/(params["radius"]*1E6) - xMargin;
        commonSetup( wg, params );
        ParaxialEquation eq;
        solver.setEquation( eq ),
        wg.setSolver( solver );
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
      default:
        clog << "Mode not recognized\n";
    }

    #ifdef VISUALIZE_PATTERN
      clog << "Visualizing waveguide intensity\n";
      visa::WindowHandler plots;
      plots.addPlot("Intensity");
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
