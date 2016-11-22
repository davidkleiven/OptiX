#include "incidentAngleSweep.hpp"
#include <string>

using namespace std;
int main( int argc, char** argv )
{
  bool useAlcohol = false;
  bool realTimeVisualize = false;
  bool saveVisualizations = false;
  for ( unsigned int i=1;i<argc;i++ )
  {
    string arg(argv[i]);
    if ( arg.find("--usealc") != string::npos )
    {
      useAlcohol = true;
    }
    else if ( arg.find("--visualize") != string::npos )
    {
      realTimeVisualize = true;
    }
    else if ( arg.find("--saveVis") != string::npos )
    {
      saveVisualizations = true;
    }
    else if ( arg.find("--help") != string::npos )
    {
      cout << "Usage: ./incidentAngleSweep.out [--useAlc --help]\n";
      cout << "help: Print this message\n";
      cout << "usealc: Use alcohol inside. Refractive indices is set for wavelength (1.57 A)\n";
      cout << "visualize: Run with real time visualization\n";
      cout << "saveVis: Save screenshot to files\n";
      return 0;
    }
    else
    {
      cout << "Unknown argument " << arg << endl;
      return 1;
    }
  }
  IncidentAngleSweep simulation;
  Visualizer vis;
  double width = 68.9;
  simulation.setWavelength( 0.1569 );
  simulation.setCladdingSilicon();
  if ( useAlcohol )
  {
    simulation.setAlcoholInside();
  }
  simulation.setWidth( width );
  simulation.setTransverseDisc( -width, 2.0*width, 1000);
  simulation.setLongitudinalDisc( 0.0, 5E6, 10000 );
  simulation.setIncAngles( -0.2, 0.2, 100 );
  simulation.setFFTSignalLength(32768);
  //simulation.saveIndx( 50 );

  if ( realTimeVisualize )
  {
    vis.setColorMax( 2.0 );
    vis.init();
    simulation.setVisualizer( vis );
    if ( saveVisualizations )
    {
      simulation.savePic("Movie");
    }
  }

  simulation.solve();


  string fname("data/angleSweep");
  simulation.save( fname );
  return 0;
}
