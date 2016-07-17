#include "dielectricSlab.h"
#include <iostream>
#include <fstream>
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

const string getArg( const std::string& arg )
{
  // Find position of =
  size_t equalPos = arg.find( "=" );
  if ( equalPos == string::npos ) return arg;
  return arg.substr( equalPos+1 );
}

/***************************** MAIN FUNCTION *********************************/
int main( int argc, char** argv )
{
  if ( argc != 4 )
  {
    cout << "Usage: ./main.out --odir=<directory of outfile> --angle=<incident angle> --mode=<steady, transient or pulse>\n";
    return 1;
  } 
  meep::initialize mpi(argc, argv);

  // ======================================================================= //
  // Parse command line arguments
  string outdir("./");
  double angle = 0.0;
  string mode("transient");
  bool odirParsed = false;
  bool angleParsed = false;
  bool modeParsed = false;
  for ( unsigned int i=1;i<argc;i++ )
  {
    stringstream ss;  
    ss << argv[i];
    if ( ss.str().find("--odir=") != string::npos )
    {
      outdir = getArg( ss.str() );
      odirParsed = true;
    }      
    else if ( ss.str().find("--angle=") != string::npos )
    {
      const string arg = getArg( ss.str() );
      ss.str("");
      ss.clear();
      ss << arg;
      ss >> angle;
      angleParsed = true;
    } 
    else if( ss.str().find("--mode=") != string::npos )
    {
      mode = getArg( ss.str() );
      modeParsed = true;
    }
    else
    {
      cout << "Unknown argument " << ss.str() << endl;
    }
    ss.str("");
    ss.clear();
  }
      
  // ======================================================================= //
  // Perform input consistency check
  bool quit = false;
  if ( !odirParsed )
  {
    cout << "Out directory was not parsed properly...\n";
    quit = true;
  }
  if ( !angleParsed )
  {
    cout << "Incident angle was not parsed properly...\n";
    quit = true;
  }
  if ( !modeParsed )
  {
    cout << "Mode was not parsed properly...\n";
    quit = true;
  }
  if ( quit ) return 1;

  if ( ( mode != "transient" ) && ( mode != "steady" ) && ( mode != "pulse" ) )
  {
    cout << "Unknown mode " << mode << endl;
    return 1;
  }

  if ( ( angle < 0.0 ) || ( angle > 90.0 ) )
  {
    cout << "Angle out of range. Has to be in the interval [0,90]\n";
    return 1;
  }

  try
  {
    string testname(outdir);
    testname += "/testfile.txt";
    ofstream of( testname.c_str() );
    string msg("Could not open file ");
    msg += testname;
    if ( !of.good() ) throw ( invalid_argument( msg ) );
  }
  catch ( invalid_argument &exc )
  {
    cout << exc.what() << endl;
    return 1;
  }
  catch (...)
  {
    cout << "An unexpected exception occured...\n";
    return 1;
  }
  // ======================================================================= //
  // Print the parameters

  cout << "Running simulation with:\n";
  cout << "Mode: " << mode << endl;
  cout << "Incident angle: " << angle << endl;
  cout << "Out directory: " << outdir << endl;

  // ======================================================================= //
  
  double resolution = 10.0;
  double freq = 0.5;

  // Initialize geometry
  DielectricSlab geometry(resolution);
  geometry.setEpsHigh(1.0);
  geometry.setEpsLow(2.25);
  double ratio = sqrt( geometry.getEpsHigh()/geometry.getEpsLow() );
  if ( ratio < 1.0 ) 
  {
    cout << "Critical angle: " << asin( ratio )*180.0/PI << endl;
  }
 
  // Compute kx
  double k = 2.0*PI*freq*sqrt( geometry.getEpsLow() );
  geometry.setKx( k*sin( angle*PI/180.0 ) );

  // Initialize source
  double width = 0.0;                   // MEEP default: 0.0
  double start_time = 0.0;              // MEEP default: 0.0
  double end_time = meep::infinity;     // MEEP default: meep::infinity
  double slowness = 10.0;               // MEEP default: 3.0
  
  meep::continuous_src_time src( freq, width, start_time, end_time, slowness  );
  meep::gaussian_src_time gaussSrc( freq, 0.4*freq );
  
  double tRelax = 4.0*geometry.getYsize(); // Time for transients to relax

  if ( mode == "transient" )
  {
    tRelax = 0.0;
  }
  double tEnd = 2.0*geometry.getYsize();   // Time to output fields

  try
  {
    geometry.addSourceVol();
    geometry.addStructure();
    geometry.addField();
    if ( mode == "pulse" )
    {
      geometry.addSource( gaussSrc, meep::Ez );
    }
    else
    {
      geometry.addSource( src, meep::Ez );  
    }
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

  geometry.getField().set_output_directory( outdir.c_str() );


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
      h5filename=outdir+"/"+h5filename+".h5";
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
