#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <cstring>
#include "readCSVdata.h"
#define DELTA_CLOSE 1E-6

using namespace std;

// ------------ MAIN FUNCTION -------------------//
int main(int argc, char** argv)
{
  if ( argc != 3 )
  {
    cout << "Usage: ./normalizeDFTFllux.cpp <infile> <bkgfile>\n";
    cout << "Arguments given:\n";
    for ( unsigned int i=0;i<argc;i++)
    {
      string arg(argv[i]);
      cout << arg << endl;
    }
    return 1;
  }

  string infile(argv[1]);
  string bkgfile(argv[2]);
 
  ReadCSVData bkg;
  ReadCSVData transRun;
  vector<double> bkgTot;
  vector<double> frequencies;
  vector<double> transmitted;
  vector<double> reflected;

  // Read the infiles
  try
  {
    bkg.read(bkgfile,3);
    transRun.read(infile, 3);
  }
  catch (runtime_error &exc)
  {
    cout << exc.what() << endl;
    return 1;
  }
  catch(...)
  {
    cout << "Unknown exception...\n";
    return 1;
  }

  // Check that the number of entris in both files are the same
  if ( bkg.numPoints() != transRun.numPoints() )
  {
    cout << "Error! Npoints bkg: "<< bkg.numPoints() <<". Npoints transmitted: " << transRun.numPoints() << endl;
    return 1;
  }

  // Normalize
  unsigned int currentBkg=0;
  bool useLastForNormalization = false;
  for ( unsigned int i=0;i<bkg.numPoints(); i++ )
  {
    frequencies.push_back(bkg.get(i,0));
    transmitted.push_back(transRun.get(i,1)/bkg.get(i,1));
    reflected.push_back(transRun.get(i,2)/bkg.get(i,2));
  }

  // Construct out filename from infile name
  size_t pos = infile.find(".");
  string ofname = infile.substr(0, pos);
  ofname += "Norm.csv";

  // Write results to file
  ofstream out(ofname.c_str());
  if ( !out.good() )
  {
    cout << "Problems when opening file " << ofname << endl;
    return 1;
  }

  out << "# Normalized total flux\n";
  out << "# Infile: " << infile << endl;
  out << "# Bkgfile: " << bkgfile << endl;
  out << "# Angle, Normalized transmitted flux, Normalized reflected\n";
  for ( unsigned int i=0;i<frequencies.size();i++)
  {
    out << frequencies[i] << "," << transmitted[i] << "," << reflected[i] << "\n";
  }
  out.close();
 
  cout << "Normalized flux written to " << ofname << endl;
  return 0;
}

