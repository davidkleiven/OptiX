
#include <iostream>
#include "meep.hpp"
#include <cmath>
#include <complex>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>
#include <stdexcept>
#include <jsoncpp/json/writer.h>

using namespace std;

/* GLOBAL PARAMETERS */
double EPS_CLAD = 1.0;

/* MAIN FUNCTION */
int main(int argc, char** argv)
{
  if ( argc != 3 )
  {
    cout << "Usage: ./slabguide.out --odir=<datadir> --epsclad=<dielectric function in cladding>\n";
    return 1;
  }

  // Parse input arguments
  string odir("");
  bool epsfound = false;
  for ( unsigned int i=1;i<argc;i++ )
  {
    string arg(argv[i]);
    if ( arg.find( "--odir=" ) != string::npos )
    {
      odir = arg.substr(7);
    }
    else if ( arg.find( "--epsclad=" ) != string::npos )
    {
      stringstream ss;
      ss << argv[i];
      ss >> EPS_CLAD;
      epsfound = true;
    }
  }
  
  // Consistency check
  if ( odir == "" )
  {
    cout << "Did not find any out directory...\n";
    return 1;
  }
  else if ( !epsfound )
  {
    cout << "Did not find epsclad argument. Using default value of 1.0...\n";
  }

  // Define geometry
