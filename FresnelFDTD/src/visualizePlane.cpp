#include "dielectricSlab.h"
#include <iostream>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include "meep.hpp"

using namespace std;
const double PI = acos(-1.0);
void convertHDF5toPng( const string& h5file, const string& pngfile, const string& epsfile )
{
  string cmd("h5topng -S3 -Zc dkbluered -m -1.0 -M 1.0 -a yarg -A ");
  cmd += epsfile;
  cmd += " -o ";
  cmd += pngfile;
  cmd += " ";
  cmd += h5file;
  FILE* pipe = popen( cmd.c_str(), "r" );
  char buffer[128];
  string result("");
  if ( !pipe )
  {
    cerr << "Could not open pipe\n";
    return;
  }
  try
  {
    while( !feof(pipe) )
    {
      if( fgets( buffer, 128, pipe ) != NULL )
      {
        result += buffer;
      }
    }
  }
  catch(...)
  {
    pclose( pipe );
    return;
  }
  pclose(pipe);
}

/***************************** MAIN FUNCTION *********************************/
int main( int argc, char** argv )
{
  if ( argc != 3 )
  {
    cout << "Usage: ./main.out <outdir> <inc angle>\n";
    return 1;
  } 

  meep::initialize mpi(argc, argv);
  stringstream ss;
  char *outdir;
  outdir = argv[1];

  ss << argv[2];
  double angle;
  ss >> angle;
  
  double resolution = 10.0;
  double freq = 0.5;

  // Initialize geometry
  DielectricSlab geometry(resolution);
  geometry.setEpsHigh(0.5);
  geometry.setEpsLow(1.0);
 
  // Compute kx
  double k = 2.0*PI*freq;
  geometry.setKx( k*sin( angle*PI/180.0 ) );

  // Initialize source
  double width = 0.0;                   // MEEP default: 0.0
  double start_time = 0.0;              // MEEP default: 0.0
  double end_time = meep::infinity;     // MEEP default: meep::infinity
  double slowness = 10.0;               // MEEP default: 3.0
  
  meep::continuous_src_time src( freq, width, start_time, end_time, slowness  );
  
  double tRelax = 4.0*geometry.getYsize(); // Time for transients to relax
  double tEnd = 2.0*geometry.getYsize();   // Time to output fields

  try
  {
    geometry.addSourceVol();
    geometry.addStructure();
    geometry.addField();
    geometry.addSource( src, meep::Ez );  
  }
  catch ( std::logic_error &exc )
  {
    cout << exc.what() << endl;
    return 1;
  }
  catch (...)
  {
    cout << "An unextepcted exception occured...\n";
    return 1;
  }

  geometry.getField().set_output_directory( outdir );


  unsigned int currentPngFile = 0;
  unsigned int nOut = 700;
  double dt = tEnd/static_cast<double>(nOut);
  double nextOutputTime = tRelax;
  while ( geometry.getField().time() < tEnd+tRelax )
  {

    geometry.getField().step();
    
    if ( geometry.getField().time() > nextOutputTime )
    {
      string h5filename("visualize");
      stringstream pngfile;
      meep::h5file* file = geometry.getField().open_h5file(h5filename.c_str());
      string odir(outdir);
      h5filename=odir+"/"+h5filename+".h5";
      pngfile << outdir << "/visualize" << currentPngFile++ << ".png";
      geometry.output_hdf5( meep::Ez, file );
      delete file;

      string epsfile(outdir);
      epsfile += "/eps-000000.00.h5";
      convertHDF5toPng( h5filename, pngfile.str(), epsfile );
      nextOutputTime += dt;
    }
  }
  return 0;
}
