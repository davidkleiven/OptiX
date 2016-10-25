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
#include <complex>
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>

using namespace std;

typedef complex<double> cdouble;

enum class Source_t {PLANE, GAUSSIAN};
int main( int argc, char **argv )
{
  double R[8] = {10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 150.0}; // In mm
  bool dumpUIDstoFile = true;
  bool useStraight = false;
  bool computeFarField = true;
  bool useCylCrd = true;
  unsigned int startRun = 0;
  unsigned int endRun = 8;
  double planeWaveAngleDeg = 0.2;
  Source_t source = Source_t::PLANE;
  /*********** PARSE COMMANDLINE ARGUMENTS ************************************/
  for ( unsigned int i=1;i<argc; i++ )
  {
    string arg(argv[i]);
    if ( arg.find("--help") != string::npos )
    {
      cout << "Usage: ./curvedWG.out [--help, --run=<run number> --straight --source=<sourcetype>]\n";
      cout << "help: Print this message\n";
      for ( unsigned int j=0;j<8;j++ )
      {
        cout << j << " Radius of curvature: " << R[j] << " mm\n";
      }
      cout << "straight: Run with a straight wave guide\n";
      cout << "source: plane or gaussian. Plane is default\n";
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
    else
    {
      cout << "Unknown argument " << arg << endl;
      return 0;
    }
  }
  /****************** END COMMANDLINE ARGUMENTS *******************************/

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
  vector<unsigned int> allUIDs;

  // Start sweep
  for ( unsigned int i=startRun;i<endRun;i++ )
  {
    clog << "Running " << i-startRun+1 << " of " << endRun-startRun << endl;
    double Rcurv = R[i]*1E6;
    double zmin = 0.0;
    double zmax = Rcurv*LzOverR;
    double xmax = width+xMarginAboveAndBelow;
    double xmin = -0.5*zmax*LzOverR-xMarginAboveAndBelow;
    double wglength = 0.9*zmax;

    if ( useStraight )
    {
      zmax = 500.0E3; // 500 um
      xmin = -xMarginAboveAndBelow;
      xmax = width + xMarginAboveAndBelow;
    }

    if ( useCylCrd )
    {
      xmin = -2.0*width;
      xmax = 2.0*width;
      // Note now z refers to the azimutal angle theta and x to the radius
      zmin = 0.0;
      double wgDistance = 500E3; // Simulate 500 um
      zmax = wgDistance/Rcurv;
    }

    double stepX = (xmax-xmin)/static_cast<double>(Nz);
    double stepZ = (zmax-zmin)/static_cast<double>(Nz);
    stepX = stepX > 1.0 ? 1.0:stepX;
    stepZ = stepZ > 100.0 ? 100.0:stepZ;

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
      else
      {
        wg = new CurvedWaveGuideFD();
      }

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
        wg->computeFarField();
        clog << "done\n";
      }
      wg->extractWGBorders();
      clog << "Exporting results...\n";
      wg->save( ctl );
      wg->saveTransmission( ctl );
      ctl.save();
      clog << "Finished exporting\n";
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
