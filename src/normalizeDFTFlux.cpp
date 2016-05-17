#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>

using namespace std;

void readFluxFile(const string &fname, vector<double> &frequencies, vector<double> &totFlux)
{
  ifstream in(fname.c_str());
  if ( !in.good() )
  {
    stringstream msg;
    msg << "Problems when opening file " << in;
    throw (msg.str());
  } 
  
  string line;
  while ( getline(in, line) && (line[0] == '#')){};
  stringstream firstLine;
  char comma;
  firstLine << line;
  double freq, newX, newY;
  firstLine >> freq >> comma >> newX >> comma >> newY;
  frequencies.push_back(freq);
  totFlux.push_back( sqrt( pow(newX, 2) + pow(newY, 2) ));
  while ( in >> freq >> comma >> newX >> comma >> newY )
  {
    frequencies.push_back(freq);
    totFlux.push_back( sqrt( pow(newX, 2) + pow(newY, 2) ));
  }
  in.close();
}

// ------------ MAIN FUNCTION -------------------//
int main(int argc, char** argv)
{
  if ( argc != 4 )
  {
    cout << "Usage: ./normalizeDFTFllux.cpp <infile> <bkgfile>\n";
    cout << "Arguments given:\n";
    for ( unsigned int i=0;i<argc;i++)
    {
      string arg(argv[i]);
      cout << arg << endl;
    }
  }

  string infile(argv[1]);
  string bkgfile(argv[2]);
 
  vector<double> bkgTot;
  vector<double> frequencies;
  vector<double> inTot;
  vector<double> inFreq;

  // Read the infiles
  try
  {
    readFluxFile(bkgfile, frequencies, bkgTot);
    readFluxFile(infile, inFreq, inTot);
  }
  catch (string &str)
  {
    cout << str << endl;
  }
  catch(...)
  {
    cout << "Unknown exception...\n";
  }
  return 0;
}

