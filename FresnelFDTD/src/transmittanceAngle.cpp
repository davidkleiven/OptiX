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
  
  unsigned int currentAngle = 0;
  vector<double> angles;
  vector<double> transmittance;
  while ( currentAngle < maxAngle )
  {
    stringstream fname;
    fname << baseFolder << currentAngle;
    ReadCSVData csvReaderNorm;
    ReadCSVData csvReaderAbs;
    try 
    {
      //string fnameNorm = fname.str() + "/transmittedFluxNorm.csv";
      string fnameNorm = fname.str() + "/ezMonitorTransFourierNorm.csv";
      csvReaderNorm.read(fnameNorm, 2);
      //string fnameAbs = fname.str() + "/transmittedFlux.csv";
      string fnameAbs = fname.str() + "/ezMonitorTransFourier.csv";
      csvReaderAbs.read(fnameAbs, 2);  
    }
    catch (runtime_error &exc)
    {
      cerr << exc.what() << endl;
      currentAngle += dAngle;
      continue;
    }
    catch(...)
    {
      cerr << "Unexpected exception occured...\n";
      currentAngle += dAngle;
      continue;
    }

    // Compute weighted sum
    double sumAbs = 0.0;
    double weightedSum = 0.0;
    if ( csvReaderNorm.numPoints() != csvReaderAbs.numPoints() )
    {
      cerr << "Different number of points in normalized and absolute files...\n";
      continue;
    }
    for ( unsigned int i=0;i<csvReaderNorm.numPoints();i++ )
    {
      sumAbs += csvReaderAbs.get(i,1);
      weightedSum += csvReaderAbs.get(i,1)*csvReaderNorm.get(i,1);
    }
    angles.push_back(currentAngle);
    transmittance.push_back(weightedSum/sumAbs); 
    currentAngle += dAngle;
  }


  // Write results to file
  string ofname("dataPlane/transmittance.csv");
  ofstream ofile(ofname.c_str());
  
  if ( !ofile.good() )
  {
    cerr << "Error when opening file " << ofname << endl;
    return 1;
  }
  ofile << "# Transmittance values\n";
  ofile << "# Angle (deg), Transmittance\n";
  for ( unsigned int i=0;i<transmittance.size();i++ )
  {
    ofile << angles[i] << "," << transmittance[i] << "\n";
  }
  ofile.close();
  return 0;
}
  
