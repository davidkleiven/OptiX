#include <iostream>
#include "cladding.hpp"
#include "curvedWaveGuide2D.hpp"
#include "crankNicholson.hpp"
#include "controlFile.hpp"
#include <complex>
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>

using namespace std;

typedef complex<double> cdouble;

int main( int argc, char **argv )
{
  double R[9] = {10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 150.0, 200.0}; // In mm
  bool dumpUIDstoFile = true;
  unsigned int startRun = 0;
  unsigned int endRun = 9;
  /*********** PARSE COMMANDLINE ARGUMENTS ************************************/
  for ( unsigned int i=1;i<argc; i++ )
  {
    string arg(argv[i]);
    if ( arg.find("--help") != string::npos )
    {
      cout << "Usage: ./curvedWG.out [--help, --run=<run number>]\n";
      cout << "help: Print this message\n";
      for ( unsigned int j=0;j<9;j++ )
      {
        cout << j << " Radius of curvature: " << R[j] << " mm\n";
      }
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
    double stepX = (xmax-xmin)/static_cast<double>(Nz);
    double stepZ = (zmax-zmin)/static_cast<double>(Nz);
    stepX = stepX > 1.0 ? 1.0:stepX;

    ControlFile ctl("data/singleCurvedWG"); // File for all parameters and settings

    try
    {
      clog << "Initializing simulation...";

      CurvedWaveGuideFD wg;
      wg.setRadiusOfCurvature( Rcurv );
      wg.setWidth( width );
      wg.setWaveLength( 0.1569 );
      wg.setCladding( cladding );
      wg.setTransverseDiscretization(xmin,xmax,stepX);
      wg.setLongitudinalDiscretization(zmin,zmax,stepZ);
      CrankNicholson solver;
      wg.setSolver(solver);
      wg.setBoundaryConditions();
      clog << " done\n";
      clog << "Solving linear system... ";
      wg.solve();
      clog << "done\n";
      clog << "Computing transmission... ";
      wg.computeTransmission( (zmax-zmin)/static_cast<double>(nPointsTransmission) );
      clog << "done\n";
      clog << "Exporting results...\n";
      wg.save( ctl );
      wg.saveTransmission( ctl );
      ctl.save();
      clog << "Finished exporting\n";

    }
    catch ( exception &exc )
    {
      cerr << exc.what() << endl;
      return 1;
    }
    catch (...)
    {
      cerr << "An unrecognized exception occured!\n";
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
