#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <cstring>
#define DELTA_CLOSE 1E-6

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

  // Normalize
  if ( bkgTot.size() != inTot.size() )
  {
    cout << "Different number of samplesin bkgfile and infile\n";
    return 1;
  }
  for ( unsigned int i=0;i<inTot.size(); i++ )
  {
    if ( abs(inFreq[i] - frequencies[i]) > DELTA_CLOSE )
    {
      cout << "Frequencies in bkgfile and infile does not match\n";
      return 1;
    }
    inTot[i] /= bkgTot[i];
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
  out << "# Frequency, Normalized total flux\n";
  for ( unsigned int i=0;i<inTot.size();i++)
  {
    out << inFreq[i] << "," << inTot[i] << "\n";
  }
  out.close();
 
  cout << "Normalized flux written to " << ofname << endl;
  return 0;
}

