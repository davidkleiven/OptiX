#include "fresnelPropagator.hpp"
#include "exitFieldSource.hpp"
#include <string>
#include <jsoncpp/json/reader.h>
#include <fstream>

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
  string paramFname("");
  /********* PARSE COMMANDLINE ARGUMENTS **************************************/
  for ( unsigned int i=1;i<argc;i++ )
  {
    string arg(argv[i]);
    if ( arg.find("--params=") != string::npos )
    {
      paramFname = arg.substr(10);
    }
    else if ( arg.find("--help") != string::npos )
    {
      cout << "Usage: ./fresnelProp.out --params=<fname> --help\n";
      cout << "phase: HDF5 file containing the phase (created by Armadillo)\n";
      cout << "amp: HDF5 file containing the amplitude (created by Armadillo)\n";
      cout << "demo: Run the demo with a step function amplitude\n";
      cout << "help: Print this message\n";
      return 0;
    }
    else
    {
      cout << "Unknown argument " << arg << endl;
      return 1;
    }
  }

  /********* END COMMANDLINE ARGUMENTS ****************************************/
  StepSource source;
  source.xmin = -30.0;
  source.xmax = 30.0;
  ExitFieldSource efSource;

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

  double xmin = params["xmin"].asDouble();
  double xmax = params["xmax"].asDouble();
  double dz = params["dz"].asDouble();
  unsigned int nTrans = params["nTransversePoints"].asInt();
  unsigned int nSteps = params["nPropagationSteps"].asInt();

  double wav = params["wavelength"].asDouble();
  bool runDemo = params["rundemo"].asBool();

  if ( !runDemo )
  {
    efSource.load( params["phase"].asString(), params["amp"].asString() );
    efSource.setDiscretization( xmin, xmax );
  }

  FresnelPropagator propagator;
  try
  {
    propagator.setWavelength(wav);
    propagator.setTransverseDiscretization( xmin, xmax, nTrans );
    if ( runDemo )
    {
      propagator.setInitialConditions( source );
    }
    else
    {
      propagator.setInitialConditions( efSource );
    }

    propagator.setStepsize( dz );
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
