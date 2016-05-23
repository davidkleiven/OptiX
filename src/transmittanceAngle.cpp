#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include "readCSVdata.h"

/**
* This program takes two command line arguments
* Arg 1 - basename for data directory
* Arg 2 - change in angle
*
* Example:
* ./transmittanceAngle.out dataPlane/Inc 10
*
* This will then read data from the folders dataPlane/Inc0, dataPlane/Inc10, dataPlane/Inc15 ...
**/

using namespace std;
int main(int argc, char** argv)
{
  if ( argc != 3 )
  {
    cout << "Usage ./transmittanceAngle.out <ddirBase> <deltaAngle>\n";
    return 1;
  }

  string baseFolder(argv[1]);
  stringstream ss;
  ss << argv[2];
  unsigned int dAngle;
  ss >> dAngle;
  const unsigned int maxAngle = 90.0;
  
  unsigned int angle = 0;
  while ( angle < maxAngle )
  {
    stringstream fname;
    fname << baseFolder << "/" << angle;
    ReadCSVData csvReaderNorm;
    ReadCSVData csvReaderAbs;
    try 
    {
      string fnameNorm = fname.str() + "transmittedFluxNorm.csv";
      csvReaderNorm.read(fnameNorm, 2);
      string fnameAbs = fname.str() + "transmittedFlux.csv";
      csvReaderAbs.read(fnameAbs, 2);  
    }
    catch (runtime_error &exc)
    {
      cerr << exc.what() << endl;
    }
    catch(...)
    {
      cerr << "Unexpected exception occured...\n";
    }

    // Compute weighted sum
    double sumAbs = 0.0;
    double weightedSum = 0.0;
  }
  return 0;
}
  
