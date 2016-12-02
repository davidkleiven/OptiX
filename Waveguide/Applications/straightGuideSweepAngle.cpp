#include "incidentAngleSweep.hpp"
#include <string>

using namespace std;
int main( int argc, char** argv )
{
  bool useAlcohol = false;
  bool realTimeVisualize = false;
  bool saveVisualizations = false;
  double delta = -1.0;
  double beta = -1.0;
  bool overrideMaterialProp = false;
  double waveguideLength = 3.0; // In mm
  bool useUID = true;
  string dfolder("data");
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
    else if ( arg.find("--deltaBeta=") != string::npos )
    {
      stringstream ss;
      ss << arg.substr(12);
      char comma;
      ss >> delta >> comma >> beta;
      overrideMaterialProp = true;
    }
    else if ( arg.find("--wglength=") != string::npos )
    {
      stringstream ss;
      ss << arg.substr(11);
      ss >> waveguideLength;
    }
    else if ( arg.find("--dfolder=") != string::npos )
    {
      dfolder = arg.substr(10);
    }
    else if ( arg.find("--noUID") != string::npos )
    {
      useUID = false;
    }
    else if ( arg.find("--help") != string::npos )
    {
      cout << "Usage: ./incidentAngleSweep.out [--useAlc --help]\n";
      cout << "help: Print this message\n";
      cout << "usealc: Use alcohol inside. Refractive indices is set for wavelength (1.57 A)\n";
      cout << "visualize: Run with real time visualization\n";
      cout << "saveVis: Save screenshot to files\n";
      cout << "deltaBeta: delta,beta. This is not recommended as it overrides the material properties!\n";
      cout << "wglength: Length of the waveguide in mm\n";
      cout << "dfolder: Data folder\n";
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
  double width = 69.8;
  double energy = 10000.0; //ev
  double lambda = 12.398*100.0/energy;
  cout << "Wavelength: " << lambda << " nm\n";
  cout << "Data folder: " << dfolder << endl;
  simulation.setWavelength( lambda );
  if ( overrideMaterialProp )
  {
    simulation.setCladdingDeltaBeta( delta, beta );
  }
  else
  {
    simulation.setCladdingSilicon( energy );
  }
  if ( useAlcohol )
  {
    simulation.setEthylenGlycolInside( energy );
    //simulation.setAlcoholInside( energy );
  }

  if ( !useUID )
  {
    simulation.turnOffUID();
  }
  simulation.setWidth( width );
  simulation.setTransverseDisc( -width, 2.0*width, 1000);
  simulation.setLongitudinalDisc( 0.0, waveguideLength*1E6, 10000 );
  clog << "Waveguide length: " << waveguideLength << " mm\n";
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

  string fname = dfolder +"/angleSweep";
  simulation.save( fname );
  return 0;
}
