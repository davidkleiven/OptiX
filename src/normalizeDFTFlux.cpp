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

void readFluxFile(const string &fname, vector<double> &angles, vector<double> &totFlux)
{
  ifstream in(fname.c_str());
  if ( !in.good() )
  {
    stringstream msg;
    msg << "Problems when opening file " << fname;
    throw (msg.str());
  } 
  
  string line;
  while ( getline(in, line) && (line[0] == '#')){};
  stringstream firstLine;
  char comma;
  firstLine << line;
  double freq, newP, angle;
  firstLine >> freq >> comma >> angle >> comma >> newP;
  angles.push_back(angle);
  totFlux.push_back( newP );
  while ( in >> freq >> comma >> angle >> comma >> newP )
  {
    angles.push_back( angle );
    totFlux.push_back( newP );
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
    return 1;
  }
  catch(...)
  {
    cout << "Unknown exception...\n";
    return 1;
  }

  // Normalize
  unsigned int currentBkg=0;
  bool useLastForNormalization = false;
  for ( unsigned int i=0;i<inTot.size(); i++ )
  {
    inTot[i] = inTot[i]/bkgTot[i];
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
  out << "# Angle, Normalized total flux\n";
  for ( unsigned int i=0;i<inTot.size();i++)
  {
    out << inFreq[i] << "," << inTot[i] << "\n";
  }
  out.close();
 
  cout << "Normalized flux written to " << ofname << endl;
  return 0;
}

