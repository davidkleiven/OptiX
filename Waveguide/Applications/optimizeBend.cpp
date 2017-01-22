#include <pei/dialogBox.hpp>
#include "bendOptimizer.hpp"
#include <map>
#include <string>
#include <stdexcept>
#include <ctime>
#include <cstdlib>

using namespace std;

int main( int argc, char** argv )
{
  srand( time(0) );
  cout << "Downsampling Z should not be modified as this will only slow down the simulation!\n";

  bool usegui = true;
  for ( unsigned int i=1;i<argc;i++ )
  {
    string arg(argv[i]);
    if ( arg.find("--nogui") != string::npos )
    {
      usegui = false;
    }
  }

  map<string,double> params;
  params["wavelength"] = 0.1569;
  params["delta"] = 4.14E-5;
  params["beta"] = 4.45E-6;
  params["stepX"] = 0.5;
  params["stepZ"] = 80.0;
  params["downSamplingZ"] = 1;
  params["downSamplingX"] = 10;
  params["width"] = 100.0;
  params["nWaveguides"] = 3;
  params["Rmin"] = 0.01;
  params["maxIter"] = 1000;
  params["totalDeflectionAngle"] = 15.0;

  if ( usegui )
  {
    pei::DialogBox box( params );
    box.show();
  }

  BendOptimizer bend;
  try
  {
    bend.init( params );
    bend.optimize();
    bend.save( "data/optimalBend.json" );
  }
  catch ( exception &exc )
  {
    cout << exc.what() << endl;
    return 1;
  }
  catch (...)
  {
    cout << "An unrecognized exception occured!\n";
    return 1;
  }
  return 0;
}
