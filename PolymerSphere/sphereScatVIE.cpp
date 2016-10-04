#include <cstdlib>
#include <iostream>
#include <libbuff.h>
#include <complex>
#include <cmath>
#include <jsoncpp/json/writer.h>
#include <fstream>
#include <string>
#include <cassert>
#include <set>
#include <stdexcept>
#include <cstdlib>
#include <ctime>

#define UID_MAX_DIGITS 6
using namespace std;

typedef complex<double> cdouble;

int main(int argc, char **argv)
{
  unsigned int UID_MAX = pow(10, UID_MAX_DIGITS);

  string help("Usage: .sphereScatVIE.out --kR=<size param> --geo=<geomatry filename>\n");
  help += "kR - size parameter: k=wave number R=radius of the sphere\n";
  help += "geo: geometry filename in accordance with the BUFF-EM documentation\n";

  double kR = -1.0;
  string geofile("");

  // Parse command line arguments
  for ( int i=1;i<argc;i++ )
  {
    string arg(argv[i]);
    if ( arg.find("--kR=") != string::npos )
    {
      stringstream kRss;
      kRss << arg.substr(5);
      kRss >> kR;
    }
    else if ( arg.find("--geo=") != string::npos )
    {
      geofile = arg.substr(6);
    }
    else if ( arg.find("--help") != string::npos )
    {
      cout << help << endl;
      return 0;
    }
    else
    {
      cerr << "Unknown argument " << arg << endl;
      return 1;
    }
  }

  // Verify that the correct arguments where given
  if ( kR < 0.0 )
  {
    cout << "No kR given.\n";
    return 0;
  }
  else if ( geofile == "" )
  {
    cout << "No geometry file given\n";
    return 0;
  }

  srand(time(0));
  unsigned int uid = rand()%UID_MAX;
  clog << "UID: " << uid << endl;
  buff::SWGGeometry geo = buff::SWGGeometry(geofile.c_str());

  // Define incident field
  double kHat[3] = {0.0,0.0,1.0};
  double omega = kR;
  cdouble E0[3];
  E0[0].real(0.0);
  E0[0].imag(0.0);
  E0[1].real(1.0);
  E0[1].imag(0.0);
  E0[2].real(0.0);
  E0[2].imag(0.0);
  PlaneWave incField(E0, kHat);
  clog << "Using plane wave source\n";

  // Define detector position
  const double detectorPosition = 1E3;
  const double deviationMax = 1.5*detectorPosition;
  const unsigned int nDetectorPixelsInEachDirection = 80;

  // Allocate memory
  HMatrix* matrix = geo.AllocateVIEMatrix();
  HVector *rhs = geo.AllocateRHSVector();

  clog << "Assembling VIE matrix...";
  geo.AssembleVIEMatrix(omega, matrix);
  clog << " done\n";

  clog << "Computing LU decomposition...";
  matrix->LUFactorize();
  clog << "done\n";

  clog << "Assembling rhs vector...";
  geo.AssembleRHSVector(omega, &incField, rhs);
  clog << "done\n";

  clog << "Solving system...";
  int info = matrix->LUSolve(rhs);
  clog << "done\n";
  delete matrix;

  // Evaluate fields quantities
  HMatrix* evaluatedFields = new HMatrix(nDetectorPixelsInEachDirection, 6, LHM_COMPLEX );
  HMatrix *Xpoints = new HMatrix(nDetectorPixelsInEachDirection, 3);
  for ( unsigned int i=0;i<nDetectorPixelsInEachDirection;i++)
  {
    double y = -deviationMax + 2.0*deviationMax*static_cast<double>(i)/static_cast<double>(nDetectorPixelsInEachDirection-1);
    Xpoints->SetEntry(i, 0, 0.0);
    Xpoints->SetEntry(i, 1, y);
    Xpoints->SetEntry(i, 2, detectorPosition);
  }

  stringstream solfname;
  Json::Value base;
  solfname << "data/solution"<<uid<<".h5";
  rhs->ExportToHDF5(solfname.str().c_str(), "SolutionVec");
  clog << "Solution vector written to " << solfname.str() << endl;

  clog << "Exporting currents currents...";
  stringstream currentFile;
  currentFile << "data/currentDistribution" << uid << ".pp";
  geo.PlotCurrentDistribution(currentFile.str().c_str(), rhs, "current");
  clog << "done\n";

  clog << "Evaluating fields...";
  stringstream fieldFile;
  fieldFile << "data/ehFields"<<uid<<".h5";
  evaluatedFields = geo.GetFields( NULL, rhs, omega, Xpoints, evaluatedFields );
  stringstream xpointFile;
  xpointFile << "data/xPointCenter" << uid << ".h5";
  base["XpointsCenter"] = xpointFile.str();
  Xpoints->ExportToHDF5(xpointFile.str().c_str(), "xPos");
  evaluatedFields->ExportToHDF5(fieldFile.str().c_str(),"Fields");
  base["FieldCenter"] = fieldFile.str();

  // Other parameters
  cdouble epsilon;
  buff::SWGVolume* vol = geo.GetObjectByLabel("sphere");
  if ( vol != NULL )
  {
    double position[3] = {0.0,0.0,0.0};
    epsilon = vol->SVT->Evaluate(omega, position);
  }

  base["Detector"]["z"] = detectorPosition;
  base["Detector"]["min"] = -deviationMax;
  base["Detector"]["max"] = deviationMax;
  base["Detector"]["pixels"] = nDetectorPixelsInEachDirection;
  base["kR"] = kR;
  base["eps"]["real"] = real(epsilon);
  base["eps"]["imag"] = imag(epsilon);
  base["nodes"] = rhs->N;
  base["UID"] = uid;
  Json::StyledWriter sw;

  stringstream overviewFile;
  overviewFile << "data/overview"<<uid<<".json";
  ofstream overview(overviewFile.str().c_str());
  if ( !overview.good() )
  {
    clog << "Could not open file " << overviewFile.str() << endl;
    return 1;
  }
  overview << sw.write(base) << endl;
  overview.close();
  clog << "Overview file written to " << overviewFile.str() << endl;

  delete rhs;
  delete evaluatedFields;
  delete Xpoints;

  if ( matrix != NULL )
  {
    // Just in case...
    delete matrix;
  }

  return 0;
}
