#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <cstring>
#include "readCSVdata.h"
#include <jsoncpp/json/reader.h>
#include <jsoncpp/json/writer.h>
#define DELTA_CLOSE 1E-6

using namespace std;

// ------------ MAIN FUNCTION -------------------//
int main(int argc, char** argv)
{
  string id("[normDFT]");
  if ( argc != 3 )
  {
    cout << id << "Usage: ./normalizeDFTFllux.cpp <infile> <bkgfile>\n";
    cout << id << "Arguments given:\n";
    for ( unsigned int i=0;i<argc;i++)
    {
      string arg(argv[i]);
      cout << arg << endl;
    }
    return 1;
  }

  string infile(argv[1]);
  string bkgfile(argv[2]);
 
  Json::Reader bkg;
  Json::Reader transRun;
  
  ifstream ibkg(bkgfile.c_str());
  if ( !ibkg.good() )
  {
    cerr << id << "Could not open file " << bkgfile << endl;
    return 1;
  }
  Json::Value rootBkg;
  bkg.parse(ibkg, rootBkg);
  ibkg.close();

  ifstream itrans(infile.c_str());
  if ( !itrans.good() )
  {
    cerr << id << "Could not open file " << infile << endl;
    return 1;
  }

  Json::Value rootTrans;
  transRun.parse( itrans, rootTrans );
  itrans.close(); 

  // Check that the number of entris in both files are the same
  if ( rootBkg["incidentAngle"].size() != rootTrans["incidentAngle"].size() )
  {
    cerr << id  << "Error! Different number of points in the background run and the transmission run\n";
    return 1;
  }

  // Normalize
  unsigned int currentBkg=0;
  bool useLastForNormalization = false;
  for ( unsigned int i=0;i<rootBkg["transmitted"].size(); i++ )
  {
    rootTrans["transmitted"][i] = rootTrans["transmitted"][i].asDouble()/rootBkg["transmitted"][i].asDouble();
    rootTrans["reflected"][i] = rootTrans["reflected"][i].asDouble()/rootBkg["reflected"][i].asDouble();
  }

  // Construct out filename from infile name
  size_t pos = infile.find(".");
  string ofname = infile.substr(0, pos);
  ofname += "Norm.json";

  rootTrans["infile"] = infile;
  rootTrans["bkgfile"] = bkgfile;
  Json::FastWriter fw;
  // Write results to file
  ofstream out(ofname.c_str());
  if ( !out.good() )
  {
    cout << id << "Problems when opening file " << ofname << endl;
    return 1;
  }

  out << fw.write(rootTrans) << endl; 
  out.close(); 
  cout << id << "Normalized flux written to " << ofname << endl;
  return 0;
}

