
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
const double YSIZE = 12.0;
const double PML = 4.0;
const double FREQ = 1.0;
const double DFREQ = 0.5*FREQ;
const double KX = 0.8;
const double XSIZE = 2.0;
const double RESOLUTION = 40.0;
unsigned int NFREQ = 5;

/* PARAMETERS DERIVED FROM THE GLOBAL */ const double CENTER = YSIZE/2.0;

/* EPSION FUNCTION */
double dielectric( const meep::vec &pos )
{
  if (( pos.y() > L+CENTER ) || ( pos.y() < CENTER-L ))
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
      ss << arg.substr(10);
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
  //const meep::symmetry sym = meep::mirror(meep::Y, vol);
  //meep::structure struc(vol, dielectric, meep::pml( PML ), sym);
  meep::structure struc(vol, dielectric, meep::pml( PML, meep::Y ));
  meep::fields field(&struc);
  field.set_output_directory( odir.c_str() );
  field.output_hdf5( meep::Dielectric, vol.surroundings() );
  meep::vec sourceLoc(2.1342, CENTER); // Random x-coordinate
  meep::gaussian_src_time source(FREQ, 0.5*FREQ);
  field.use_bloch( meep::X, KX );
  field.add_point_source( meep::Ez, source, sourceLoc ); 

  // Compute time required to propagate over the entire domain
  double tProp = XSIZE; // Speed of light c=1
  
  // Define output quantities
  Json::Value monitorInside(Json::arrayValue);
  Json::Value monitorsInsidePos(Json::arrayValue);
  Json::Value monitorOutside(Json::arrayValue);
  Json::Value monitorsOutsidePos(Json::arrayValue);
  double monitorX = XSIZE -1.0;
  unsigned int nMonitors = 11;
  double dyInside = 2.0*static_cast<double>(L)/static_cast<double>(nMonitors-1); 
  double dyOutside = static_cast<double>(L)/static_cast<double>(nMonitors);
  // Initialize arrays
  for ( unsigned int i=0;i<nMonitors;i++)
  {
    monitorInside.append(0.0);
    monitorOutside.append(0.0);
    monitorsInsidePos.append( CENTER-L+static_cast<double>(i)*dyInside);
    monitorsOutsidePos.append( CENTER + L + static_cast<double>(i)*dyOutside);
  }
  
  // Setup DFT point to collect the frequency
  vector< complex<double> > fieldMonitor;
  while ( field.time() < source.last_time()+10.0*tProp )
  {
    field.step();
    meep::vec fieldPos(XSIZE/2.2, CENTER);
    fieldMonitor.push_back( field.get_field( meep::Ez, fieldPos) );
  }
  
  // Use harminv
  vector< complex<double> > amps(NFREQ);
  vector<double> freq_re(NFREQ);
  vector<double> freq_im(NFREQ);
  int status = meep::do_harminv( &fieldMonitor[0], fieldMonitor.size(), field.dt, FREQ-DFREQ, FREQ+DFREQ, NFREQ, &amps[0], &freq_re[0], &freq_im[0] );
  // Collect the frequency
  Json::Value freqs(Json::arrayValue);
  Json::Value qFactor(Json::arrayValue);
  for ( unsigned int i=0;i<freq_re.size();i++ )
  {
    freqs.append( freq_re[i] );
    qFactor.append( -freq_re[i]/(2.0*freq_im[i] ));
  }
  

  // Get the field
  for ( unsigned int i=0;i<nMonitors;i++ )
  { 
    complex<double> Ez = field.get_field( meep::Ez, meep::vec(monitorX, monitorsInsidePos[i].asDouble()) );
    monitorInside[i] = abs(Ez);
    Ez = field.get_field( meep::Ez, meep::vec(monitorX, monitorsOutsidePos[i].asDouble()) );
    monitorOutside[i] = abs(Ez);
  }

  Json::Value base;
  base["geometry"]["EpsClad"] = EPS_CLAD;
  base["geometry"]["guideWidth"] = 2.0*L;
  base["monitor"]["inside"]["data"] = monitorInside;
  base["monitor"]["inside"]["pos"] = monitorsInsidePos;
  base["monitor"]["outside"]["data"] = monitorOutside;
  base["monitor"]["outside"]["pos"] = monitorsOutsidePos;
  base["monitor"]["freq"] = freqs;
  base["monitor"]["qFactor"] = qFactor;
  Json::FastWriter fw;
  string ofname = odir+"/monitors.json";
  ofstream os(ofname.c_str());
  if ( !os.good() )
  {
    cout << "Problem when opening file " << ofname << endl;
    return 1;
  }
  
  os << fw.write( base ) << endl;
  os.close();
  cout << "File written to " << ofname << endl;
  // Output field in the end
  field.output_hdf5( meep::Ez, vol.surroundings() ); 
  return 0;
}
