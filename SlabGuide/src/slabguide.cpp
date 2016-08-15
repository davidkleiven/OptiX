
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
const double L = 1.0;
const double YSIZE = 5.0;
const double PML = 2.0;
const double FREQ = 1.0;
const double XSIZE = 10.0;
const double RESOLUTION = 10.0;

/* EPSION FUNCTION */
double dielectric( const meep::vec &pos )
{
  if ( pos.y() > L )
  {
    return EPS_CLAD;
  }
  return 1.0;
} 

/* MAIN FUNCTION */
int main(int argc, char** argv)
{
  meep::initialize mpi(argc, argv); // Standard beginning recommended by the MEEP authors
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
  meep::grid_volume vol = meep::vol2d( XSIZE, YSIZE, RESOLUTION );
  const meep::symmetry sym = meep::mirror(meep::Y, vol);
  meep::structure struc(vol, dielectric, meep::pml( PML ), sym);
  meep::fields field(&struc);
  field.set_output_directory( odir.c_str() );
  meep::vec sourceLoc(PML, 0.0);
  meep::continuous_src_time source(FREQ);
  field.add_point_source( meep::Ez, source, sourceLoc ); 

  // Compute time required to propagate over the entire domain
  double tProp = XSIZE; // Speed of light c=1
  
  while ( field.time() < 1.5*tProp )
  {
    field.step();
  }
  
  // Output field in the end
  field.output_hdf5( meep::Ez, vol.surroundings() ); 
  return 0;
}
