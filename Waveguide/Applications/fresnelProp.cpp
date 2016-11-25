#include "fresnelPropagator.hpp"
#include "exitFieldSource.hpp"
#include <string>
#include <jsoncpp/json/reader.h>
#include <fstream>
#include <H5Cpp.h>
#include <cmath>

class StepSource
{
public:
  double operator()( double x ) const
  {
    if ( ( x > xmin) && ( x < xmax) )
    {
      return 1.0;
    }
    return 0.0;
  }
  double xmin, xmax;
};

using namespace std;
int main( int argc, char** argv )
{
  const double PI = acos(-1.0);
  string paramFname("");
  double zmax=1.0;
  unsigned int nSteps=0;
  bool runDemo = false;
  double xmin, xmax;
  unsigned int Nx=0;
  /********* PARSE COMMANDLINE ARGUMENTS **************************************/
  for ( unsigned int i=1;i<argc;i++ )
  {
    string arg(argv[i]);
    if ( arg.find("--params=") != string::npos )
    {
      paramFname = arg.substr(9);
    }
    else if ( arg.find("--help") != string::npos )
    {
      cout << "Usage: ./fresnelProp.out --params=<fname> --help\n";
      cout << "params: The far field file containing both exit field amplitude and phase\n";
      cout << "zmax: Maximum position in um\n";
      cout << "nsteps: Number of steps\n";
      cout << "xdisc: xmin,xmax,Nx\n";
      cout << "demo: Run the demo with a step function amplitude\n";
      cout << "help: Print this message\n";
      return 0;
    }
    else if ( arg.find("--nsteps=") != string::npos )
    {
      stringstream ss;
      ss << arg.substr(9);
      ss >> nSteps;
    }
    else if ( arg.find("--zmax=") != string::npos )
    {
      stringstream ss;
      ss << arg.substr(7);
      ss >> zmax;
    }
    else if ( arg.find("--xdisc=") != string::npos )
    {
      stringstream ss;
      ss << arg.substr(8);
      char comma;
      ss >> xmin >> comma >> xmax >> comma >> Nx;
    }
    else
    {
      cout << "Unknown argument " << arg << endl;
      return 1;
    }
  }

  if ( nSteps == 0)
  {
    cout << "Number of steps specified is zero!\n";
    return 1;
  }

  if ( Nx == 0 )
  {
    cout << "No transverse discretization given!\n";
    return 1;
  }

  /********* END COMMANDLINE ARGUMENTS ****************************************/
  StepSource source;
  source.xmin = -30.0;
  source.xmax = 30.0;
  ExitFieldSource efSource;

  /*
  Json::Value params;
  Json::Reader reader;
  ifstream infile(paramFname.c_str());
  if ( !infile.good() )
  {
    cout << "Error when reading the parameter file\n";
    return 1;
  }
  reader.parse(infile, params);
  infile.close();
  */

  // Open HDf5 file
  /*
  double xmin = params["xmin"].asDouble();
  double xmax = params["xmax"].asDouble();
  double dz = params["dz"].asDouble();
  unsigned int nTrans = params["nTransversePoints"].asInt();
  unsigned int nSteps = params["nPropagationSteps"].asInt();

  double wav = params["wavelength"].asDouble();
  bool runDemo = params["rundemo"].asBool();
  */

  double wav;
  try
  {
    if ( !runDemo )
    {
      H5::H5File file(paramFname, H5F_ACC_RDONLY );
      H5::DataSet farField = file.openDataSet("farField");
      H5::Attribute attr = farField.openAttribute("wavenumber");
      H5::DataType type = attr.getDataType();
      attr.read(type, &wav);
      wav = 2.0*PI/wav;
      efSource.load( file );
    }
  }
  catch( H5::Exception &exc )
  {
    cout << exc.getCDetailMsg() << endl;
    return 1;
  }

  FresnelPropagator propagator;
  try
  {
    propagator.setWavelength(wav);
    if ( runDemo )
    {
      propagator.setInitialConditions( source );
    }
    else
    {
      propagator.setTransverseDiscretization( xmin, xmax, Nx );
      propagator.setInitialConditions( efSource );
    }

    propagator.setStepsize( zmax/nSteps );
    propagator.propagate( nSteps );
    string fname("data/fresnelPropTest.h5");
    propagator.save( fname );
  }
  catch ( exception &exc )
  {
    cout << exc.what() << endl;
    return 1;
  }
  catch (...)
  {
    cout << "Unexpected exception!\n";
    return 1;
  }
  return 0;
}
