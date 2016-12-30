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

int main( int argc, char **argv )
{
  double R[8] = {10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 150.0}; // In mm
  double widths[8] = {40.0, 60.0, 80.0, 100.0, 200.0, 400.0, 800.0, 1200.0}; // In nm
  double logDeltas[8] = {-6.0, -5.5, -5.25, -5.0, -4.75, -4.5, -4.25, -4.0}; // Used for sweep

  bool dumpUIDstoFile = true;
  bool useStraight = false;
  bool computeFarField = true;
  bool useCylCrd = false;
  bool useBorderTracker = false;
  bool performWidthSweep = false;
  bool deltaSweep = false;
  const unsigned int defaultWidthIndex = 6;
  unsigned int widthStart = defaultWidthIndex;
  unsigned int widthEnd = widthStart+1;
  unsigned int startRun = 0;
  unsigned int endRun = 8;
  unsigned int deltaStart = 0;
  unsigned int deltaEnd = 1;
  unsigned int Nx = 2000;
  unsigned int Nz = 6000;
  double planeWaveAngleDeg = 0.0;//0.2;
  double linearRampWidthFraction = 5.0; // Only relevant of profile = LINEAR_RAMP
  Source_t source = Source_t::PLANE;
  WGProfile_t profile = WGProfile_t::STEP;

  /*********** PARSE COMMANDLINE ARGUMENTS ************************************/
  for ( unsigned int i=1;i<argc; i++ )
  {
    string arg(argv[i]);
    if ( arg.find("--help") != string::npos )
    {
      cout << "Usage: ./curvedWG.out [--help, --run=<run number> --straight --source=<sourcetype> --disc=Nx,Nz]\n";
      cout << "help: Print this message\n";
      for ( unsigned int j=0;j<8;j++ )
      {
        cout << j << " Radius of curvature: " << R[j] << " mm\n";
      }
      cout << "straight: Run with a straight wave guide\n";
      cout << "source: plane or gaussian. Plane is default\n";
      cout << "cyl: use cylindrical coordiantes\n";
      cout << "brdtrack: Use the border tracker. Only have effect if using cartesian coordinates\n";
      cout << "widthsweep: Perform a sweep over different widths\n";
      cout << "deltasweep: Perform a sweep over different delta\n";
      cout << "disc: Number of discretization points in each direction";
      return 0;
    }
    else if ( arg.find("--run=") != string::npos )
    {
      stringstream ss;
      ss << arg.substr(6);
      ss >> startRun;
      endRun = startRun + 1;
      dumpUIDstoFile = false;
    }
    else if ( arg.find("--straight") != string::npos )
    {
      useStraight = true;
      endRun = 1;
      dumpUIDstoFile = false;
    }
    else if ( arg.find("--source=") != string::npos )
    {
      string sourceStr(arg.substr(9));
      if ( sourceStr == "gaussian" )
      {
        source = Source_t::GAUSSIAN;
      }
      else if ( sourceStr == "plane" )
      {
        source = Source_t::PLANE;
      }
      else
      {
        cout << "Unknown source type " << sourceStr << endl;
        return 0;
      }
    }
    else if ( arg.find("--cyl") != string::npos )
    {
      useCylCrd = true;
    }
    else if ( arg.find("--brdtrack") != string::npos )
    {
      useBorderTracker = true;
    }
    else if ( arg.find("--widthsweep") != string::npos )
    {
      performWidthSweep = true;
    }
    else if ( arg.find("--deltaSweep") != string::npos )
    {
      deltaSweep = true;
    }
    else if ( arg.find("--disc=") != string::npos )
    {
      stringstream ss;
      ss << arg.substr(7);
      char comma;
      ss >> Nx >> comma >> Nz;
    }
    else
    {
      cout << "Unknown argument " << arg << endl;
      return 0;
    }
  }
  /****************** END COMMANDLINE ARGUMENTS *******************************/

  // Parameters for running a sweep over radii of curvature
  double LzOverR = 0.01; // max(z)/R << 1 is a requirement
  double xMarginAboveAndBelow = 0.01E3; // In nanometers = 0.5 um

  unsigned int nPointsTransmission = 200;

  Cladding cladding;
  double delta = 4.14E-5; // Default value
  double beta = 3.45E-6;
  cladding.setRefractiveIndex(delta, beta);
  double width = 100.0; // Width of the waveguide in nm
  vector<unsigned int> allUIDs;
  bool allowUseOfBorderTracker = false; // Internally handled, do not change this

  if ( performWidthSweep )
  {
    widthStart = 0;
    widthEnd = 8;
  }

  if ( deltaSweep )
  {
    deltaStart = 0;
    deltaEnd = 8;
  }

  // Start sweep
  for ( unsigned int i=startRun;i<endRun;i++ )
  {
    clog << "Running " << i-startRun+1 << " of " << endRun-startRun << endl;
    for ( unsigned int iw=widthStart;iw<widthEnd;iw++ )
    {
      for ( unsigned int idelta=deltaStart;idelta<deltaEnd;idelta++ )
      {
        if ( deltaSweep )
        {
          delta = pow( 10.0, logDeltas[idelta] );
          cladding.setRefractiveIndex(delta, beta);
          clog << "Running width delta = " << delta << endl;
        }
        double Rcurv = R[i]*1E6;
        double zmin = 0.0;
        double zmax = Rcurv*LzOverR;
        double xmax = width+xMarginAboveAndBelow;
        double xmin = -0.5*zmax*LzOverR-xMarginAboveAndBelow;
        double wglength = 1.5*zmax;

        if ( performWidthSweep )
        {
          width = widths[iw];
        }
        clog << "Running width waveguide width: " << width << "nm\n";

        if ( useStraight )
        {
          zmax = 400.0E3; // 500 um
          xmin = -width;
          xmax = width + width;
          wglength = 1.5*zmax;
        }

        if ( useCylCrd )
        {
          xmin = -width;
          xmax = width+width;
          // Note now z refers to the azimutal angle theta and x to the radius
          zmin = 0.0;
          double wgDistance = 400E3; // Simulate 400 um

          // Hack for generating the geometrical optics plot
          //width = 1200.0;
          //xmin = -0.2*width;
          //xmax = width + 0.2*width;
          //wgDistance = 1200E3;

          zmax = wgDistance/Rcurv;
          wglength = 1.1*wgDistance;
        }

        //stepX = stepX > 1.0 ? 1.0:stepX;
        //stepZ = stepZ > 100.0 ? 100.0:stepZ;

        ControlFile ctl("data/singleCurvedWG"); // File for all parameters and settings

        CurvedWaveGuideFD *wg = NULL;
        ParaxialSource *src = NULL;
        ParaxialEquation *eq = NULL; // Cartesian coordinates
        try
        {
          clog << "Initializing simulation...";
          if ( useStraight )
          {
            wg = new StraightWG2D();
          }
          else if ( useCylCrd )
          {
            wg = new CurvedWGCylCrd();
          }
          else if ( profile == WGProfile_t::GAUSSIAN )
          {
            wg = new GaussianWG();
          }
          else if ( profile == WGProfile_t::LINEAR_RAMP )
          {
            LinearRampWG* linWG = new LinearRampWG();
            linWG->setWidthFraction( linearRampWidthFraction );
            wg = linWG;
          }
          else
          {
            wg = new CurvedWaveGuideFD();
            allowUseOfBorderTracker = true;
          }

          if ( useBorderTracker && allowUseOfBorderTracker )
          {
            xmin = -1.0*width;
            xmax = 2.0*width;
          }

          double stepX = (xmax-xmin)/static_cast<double>(Nx);
          double stepZ = (zmax-zmin)/static_cast<double>(Nz);

          double wavelength = 0.1569;
          switch ( source )
          {
            case Source_t::PLANE:
            {
              clog << "Using plane wave source\n";
              PlaneWave *psrc = new PlaneWave();
              psrc->setAngleDeg( planeWaveAngleDeg );
              src = psrc;
              break;
            }
            case Source_t::GAUSSIAN:
              clog << "Using gaussian source\n";
              GaussianBeam *gsrc = new GaussianBeam();
              gsrc->setWaist( 1.0 ); // Waist is 10 nm with a wavelength of 0.1569 nm this gives beamdivergence of 0.3 deg
              gsrc->setOrigin( width/2.0, -1E4 ); // Origin at -100.0 nm
              src = gsrc;
              break;
          }

          src->setWavelength( wavelength );

          wg->setRadiusOfCurvature( Rcurv );
          wg->setWaveguideLength( wglength );
          wg->setWidth( width );
          wg->setWaveLength( wavelength );
          wg->setCladding( cladding );
          wg->setTransverseDiscretization(xmin,xmax,stepX);
          wg->setLongitudinalDiscretization(zmin,zmax,stepZ);
          CrankNicholson solver;
          if ( useCylCrd )
          {
            CylindricalParaxialEquation *ceq = new CylindricalParaxialEquation();
            ceq->setRadiusOfCurvature( Rcurv );
            eq = ceq;
          }
          else
          {
            eq = new ParaxialEquation();
          }

          solver.setEquation( *eq );
          wg->setSolver(solver);
          wg->setBoundaryConditions( *src );
          if ( useBorderTracker )
          {
            wg->useBorderTracker();
          }
          clog << " done\n";
          clog << "Solving linear system... ";
          wg->solve();
          clog << "done\n";

          clog << "Computing transmission... ";
          wg->computeTransmission( (zmax-zmin)/static_cast<double>(nPointsTransmission) );
          clog << "done\n";

          if ( computeFarField )
          {
            clog << "Computing far fields... ";
            wg->setFarFieldAngleRange( -0.5, 0.5 );
            wg->computeFarField( 65536 );
            clog << "done\n";
          }

          if ( !useBorderTracker )
          {
            wg->extractWGBorders();
          }

          clog << "Exporting results...\n";
          #ifdef VISUALIZE_PATTERN
            //wg->saveContour( false );
          #endif
          wg->save( ctl );
          wg->saveTransmission( ctl );
          ctl.save();
          clog << "Finished exporting\n";

          #ifdef VISUALIZE_PATTERN
            clog << "Visualizing waveguide intensity\n";
            visa::WindowHandler plots;
            plots.addPlot("Intensity");
            arma::mat intensity = arma::abs( wg->getSolver().getSolution() );
            plots.get("Intensity").fillVertexArray( intensity );
            plots.show();
            for ( unsigned int i=0;i<KEEP_PLOT_FOR_SEC;i++ )
            {
              plots.show();
              clog << "Closes in " << KEEP_PLOT_FOR_SEC-i <<  "seconds \r";
              this_thread::sleep_for( chrono::seconds(1) );
            }
          #endif

          delete wg;
          delete src;
          delete eq;
        }
        catch ( exception &exc )
        {
          cerr << exc.what() << endl;
          delete wg;
          delete src;
          delete eq;
          return 1;
        }
        catch (...)
        {
          cerr << "An unrecognized exception occured!\n";
          delete src;
          delete wg;
          delete eq;
          return 1;
        }

        clog << "The simulation ended successfully\n";
        allUIDs.push_back( ctl.getUID() );
      }
    }
  }

  if ( dumpUIDstoFile )
  {
    string uidFname = "data/allUIDs.txt";
    ofstream out( uidFname.c_str() );
    for ( unsigned int i=0; i<allUIDs.size(); i++ )
    {
      out << allUIDs[i] << endl;
    }
    out.close();
    clog << "UIDs for sweep written to " << uidFname << endl;
  }
  return 0;
}
