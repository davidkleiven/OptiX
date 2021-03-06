#include <visa/visa.hpp>
#include <pei/dialogBox.hpp>
#include "multipleCurvedWG.hpp"
#include <map>
#include <string>
#include <chrono>
#include <thread>
#include <iostream>
#define KEEP_PLOT_FOR_SEC 6

using namespace std;

typedef visa::Colormaps::Colormap_t cmap_t;
int main ( int argc, char** argv )
{
  string geofile("");
  for ( unsigned int i=1;i<argc;i++ )
  {
    string arg(argv[i]);
    if ( arg.find("--geo=") != string::npos )
    {
      geofile = arg.substr( 6 );
    }
    else if ( arg.find("--help") != string::npos )
    {
      cout << "Usage: ./multiWG.out --geo=<geofile> [--help]\n";
      cout << "geo: JSON file containing the waveguide geometry\n";
      cout << "help: Print this message\n";
      return 0;
    }
    else
    {
      cout << "Unknown argument " << arg << endl;
      return 1;
    }
  }

  if ( geofile == "" )
  {
    cout << "No geometry file given!\n";
    return 1;
  }

  map<string,double> params;
  params["wavelength"] = 0.1569;
  params["delta"] = 4.14E-5;
  params["beta"] = 4.45E-6;
  params["stepX"] = 0.5;
  params["stepZ"] = 80.0;
  params["downSamplingZ"] = 10.0;
  params["downSamplingX"] = 10;
  params["width"] = 100.0;

  {
    pei::DialogBox dialog( params );
    dialog.show();
  }

  MultipleCurvedWG simulation;
  try
  {
    simulation.loadWaveguides( geofile );
    simulation.init( params );
    simulation.solve();

    simulation.save( "data/multipleWG.h5" );
  }
  catch( exception &exc )
  {
    cout << exc.what() << endl;
    return 1;
  }
  catch(...)
  {
    cout << "Unrecogniezed exception!\n";
    return 1;
  }

}
