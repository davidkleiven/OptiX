#include "fresnelPropagator.hpp"
#include "exitFieldSource.hpp"
#include <string>
#include <json/reader.h>
#include <fstream>
#include <H5Cpp.h>
#include <cmath>
#include <ctime>
#include <visa/visa.hpp>
#include <pei/dialogBox.hpp>
#include <map>
#include <thread>
#include <chrono>
#define VISUALIZE_INTENSITY
#define KEEP_PLOT_FOR_SEC 5

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
typedef visa::Colormaps::Colormap_t cmap_t;
int main( int argc, char** argv )
{
  srand(time(0));
  const double PI = acos(-1.0);
  string paramFname("");
  map<string,double> params;
  params["zmax"] = 1.0;
  params["nSteps"] = 100;
  params["xmin"] = -100.0;
  params["xmax"] = 100.0;
  params["Nx"] = 100;

  pei::DialogBox box( params );
  box.show();

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
      cout << "params: The HDF5 file containing both exit field amplitude and phase\n";
      cout << "help: Print this message\n";
      return 0;
    }
    else
    {
      cout << "Unknown argument " << arg << endl;
      return 1;
    }
  }

  if ( static_cast<int>(params["nSteps"]+0.5) == 0)
  {
    cout << "Number of steps specified is zero!\n";
    return 1;
  }

  if (  static_cast<int>( params["Nx"] + 0.5 ) == 0 )
  {
    cout << "No transverse discretization given!\n";
    return 1;
  }

  StepSource source;
  source.xmin = -30.0;
  source.xmax = 30.0;
  ExitFieldSource efSource;

  double wav;
  try
  {
      H5::H5File file(paramFname, H5F_ACC_RDONLY );
      H5::DataSet farField = file.openDataSet("farField");
      H5::Attribute attr = farField.openAttribute("wavenumber");
      H5::DataType type = attr.getDataType();
      attr.read(type, &wav);
      wav = 2.0*PI/wav;
      efSource.load( file );
  }
  catch( H5::Exception &exc )
  {
    cout << exc.getCDetailMsg() << endl;
    return 1;
  }

  FresnelPropagator propagator;
  try
  {
    propagator.setWavelength( wav );
    propagator.setTransverseDiscretization( params.at("xmin"), params.at("xmax"), params.at("Nx") );
    propagator.setInitialConditions( efSource );
    propagator.setStepsize( params.at("zmax")/params.at("nSteps") );
    propagator.propagate( params.at("nSteps") );
    string fname("data/fresnelProp");
    propagator.save( fname );

    #ifdef VISUALIZE_INTENSITY
      arma::mat intensity = propagator.getIntensity();
      //cmap_t cmap = cmap_t::VIRIDIS;
      cmap_t cmap = cmap_t::NIPY_SPECTRAL;
      visa::WindowHandler plots;
      plots.addPlot("Intensity");
      plots.get("Intensity").setCmap( cmap );
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
