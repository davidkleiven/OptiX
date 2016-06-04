#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "readCSVdata.h"
#include <stdexcept>
#define CLOSE_DELTA 1E-8

using namespace std;

typedef vector<double> dvec;
int main(int argc, char** argv)
{
  if ( argc != 3 )
  {
    cout << "Usage: ./subtractBackground.cpp <ifile> <bkgfile>\n";
    return 1;
  }

  string ifile(argv[1]);
  string bkgfile(argv[2]);

  dvec timeIn;
  dvec timeBkg;
  dvec realIn;
  dvec realBkg;
  dvec imagIn;
  dvec imagBkg;

  ReadCSVData bkgReader;
  ReadCSVData ifileReader;
  try
  {
    bkgReader.read(bkgfile, 3);
    ifileReader.read(ifile, 3);
  }
  catch ( runtime_error &exc )
  {
    cerr << exc.what() << endl;
    return 1;
  }
  catch (...)
  {
    cerr << "Unrecognized exception occured...\n";
    return 1;
  }


  // Subtract off the background
  if ( bkgReader.numPoints() != bkgReader.numPoints() )
  {
    cout << "The background time and the bkg time must have exactly the same timepoints!\n";
    return 1;
  }

  vector<double> reflectedSubtracted;
  for ( unsigned int i=0;i<bkgReader.numPoints();i++)
  {
    if ( abs( bkgReader.get(i,0) - ifileReader.get(i,0) ) > CLOSE_DELTA )
    {
      cout << "Timepoints on step " << i << "does not match!\n";
      cout << "Input file: " << ifileReader.get(i,0) << endl;
      cout << "Bkg file: " << bkgReader.get(i,0) << endl << endl;
      return 1;
    }
    reflectedSubtracted.push_back( ifileReader.get(i,2) - bkgReader.get(i,2) );
  }
  
  // Write the results to file
  ofstream os(ifile.c_str());
  if ( !os.good() )
  {
    cout << "Problems when opening file " << ifile << " for writing...\n";
    return 1;
  }

  os << "# Field after subtracting of background\n";
  os << "# Infile: " << ifile << "\n";
  os << "# Background file " << bkgfile << "\n";
  os << "# Time, Real field transmitted, Real field reflected (incident is subtracted)\n";
  for ( unsigned int i=0;i<bkgReader.numPoints();i++)
  {
    os << bkgReader.get(i,0) << "," << ifileReader.get(i,1) << "," << reflectedSubtracted[i] << "\n";
  }
  os.close();
  return 0;
}
