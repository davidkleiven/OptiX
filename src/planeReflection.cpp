#include <iostream>
#include "meep.hpp"
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>
#include <stdexcept>
#include "sincSrc.h"
#define ODIRLEN 60
//#define OUTPUT_HDF5

/**
* This program takes 3 command line arguments
* Arg 1 - out directory. Directory where the datafiles will be stored
* Arg 2 - Relative dielectriv permittivity in scattering region. The wave starts in a region with epsilon=1.
* Arg 3 - angle of incidence in degrees
* Arg 4 - Polarization (s or p)
*
* Example:
*   ./planeReflection data/planewave 2.25 30 s
*/
using namespace std;

const double EPS_LOW = 1.0;
double EPS_HIGH = 1.0;

const double XSIZE = 5.0;
const double YSIZE = 18.0;
const double PML_THICK = 6.0;
const double SOURCE_Y = YSIZE-PML_THICK - 1.0;
const double YC_PLANE = YSIZE/2.0;
//const double YC_PLANE = PML_THICK+1.0;
const double PI = acos(-1.0);
const complex<double> IMAG_UNIT(0,1.0);
const unsigned int NSTEPS = 20;
double KX = 1.0;

double dielectric(const meep::vec &pos)
{
  if ( pos.y() < YC_PLANE )
  {
    return EPS_HIGH;
  }
  return EPS_LOW;
}

complex<double> amplitude(const meep::vec &pos)
{
  return exp(IMAG_UNIT*KX*pos.x());
}

//----------------- MAIN FUNCTION ----------------------------//
int main(int argc, char **argv)
{
  cout << "This program simulates scattering of a plane wave onto a smooth surface\n";
  meep::initialize mpi(argc, argv);

  // Read command line arguments
  if ( argc != 7 )
  {
    cout << "Usage: ./planeReflection.out <out directory> <epsInScattered> <incidend angle> <polarization> <relBandwidth>";
    cout << "<number of frequencies>\n";
    cout << "The following arguments were given:\n";
    for ( unsigned int i=0;i<argc;i++ )
    {
      string arg(argv[i]);
      cout << arg << endl;
    }
    return 1;
  }

  const char* OUTDIR = argv[1];
  stringstream ss;
  ss << argv[2];
  ss >> EPS_HIGH;
  ss.clear();
  ss << argv[3];
  double angle;
  ss >> angle;
  char polarization = argv[4][0];
  double relFwidth;
  ss.clear();
  ss << argv[5];
  ss >> relFwidth;

  unsigned int nfreq;
  ss.clear();
  ss << argv[6];
  ss >> nfreq;


  // Check that angle is within range
  const double maxAngle = 90.0;
  if ( angle > maxAngle )
  {
    cout << "The incident angle is too large\n";
    cout << "Maximum angle is " << maxAngle << endl;
    return 1;
  }
  else if ( angle < 0.0 )
  {
    cout << "Negative angle given. Has to be in range [0,MAX_ANGLE)\n";
    return 1;
  }

  // Check that polarization holds a valid value
  if ( (polarization != 's') && (polarization != 'p') )
  {
    cout << "Invalid polarization value. Has to be either s or p\n";
    cout << "Value given: " << polarization << endl;
    return 1;
  }

  double freq = 0.3;
  double fwidth = freq*relFwidth;

  // Compute kx
  double k = 2.0*PI*freq;
  KX = k*sin( angle*PI/180.0 );

  const double minHeight = 2.0*PML_THICK + 4.0;

  // Verify that the size of the domain is big enough (for debugging only)
  assert ( YSIZE > minHeight );
  
  double resolution = 5.0;

  // Initialize computational cell
  meep::grid_volume vol = meep::vol2d(XSIZE, YSIZE, resolution);

  // Add line source to get plane wave
  meep::vec srcCorner1(0.0, SOURCE_Y);
  meep::vec srcCorner2(XSIZE, SOURCE_Y);
  meep::volume srcvol(srcCorner1, srcCorner2);

  // Initalize structure. Add PML in y-direction
  meep::structure srct(vol, dielectric, meep::pml( PML_THICK, meep::Y ) );

  //srct.Courant = 0.1;
  meep::fields field(&srct);
  field.use_bloch( meep::X, KX/(2.0*PI) ); // MEEP leaves out the factor 2pi (k = 1/lambda)
  field.set_output_directory(OUTDIR); 
  
  // Write dielectric function to file
  field.output_hdf5(meep::Dielectric, vol.surroundings()); 

  // Set source type. Use Gaussian, had some problems with the continous
  meep::gaussian_src_time src(freq, fwidth);
  
  meep::component fieldComp;
  if ( polarization == 's' )
  {
    fieldComp = meep::Ez;
  }
  else
  {
    fieldComp = meep::Hz;
  }

  field.add_volume_source(fieldComp, src, srcvol, amplitude);
  // Add DFT fluxplane
  const double fluxPlanePosY = 0.5*(YC_PLANE + PML_THICK);
  const double fluxRefPlanePosY = 0.5*(YC_PLANE + SOURCE_Y);
  meep::volume dftVol = meep::volume(meep::vec(0.0,fluxPlanePosY), meep::vec(XSIZE-1.0,fluxPlanePosY));
  meep::volume dftVolR = meep::volume(meep::vec(0.0,fluxRefPlanePosY), meep::vec(XSIZE-1.0,fluxRefPlanePosY)); 
  meep::dft_flux transFluxY = field.add_dft_flux_plane(dftVol, freq-fwidth/2.0, freq+fwidth/2.0, nfreq);
  meep::dft_flux fluxYReflected = field.add_dft_flux_plane(dftVolR, freq-fwidth/2.0, freq+fwidth/2.0, nfreq); 

  vector<double> fieldTransmitted; // Container for the real field component
  vector<double> timepoints; // Container for the timepoints
  vector<double> fieldReflection; // Container for the real field component at the center of the flux plane

  // Time required to propagate over the domain with the slowest speed
  double speed = 1.0/sqrt(EPS_HIGH);
  double tPropagate = field.last_source_time() + 1.2*SOURCE_Y/speed;

  // Main loop.
  double transYWidth = abs( dftVol.get_max_corner().x() - dftVol.get_min_corner().x() );
  double timeToRegisterFourier = nfreq/fwidth;

  double tEnd = timeToRegisterFourier > tPropagate ? timeToRegisterFourier:tPropagate;

  unsigned int nOut = 80; // Number of output files
  double dt = tEnd/nOut;
  double nextOutputTime = 0.0;
  string fluxXFname("fluxYReflected");

  string outdir(OUTDIR);
  if ( outdir.find("bkg") == string::npos )
  {
    // Did not find bkg in the directory name. Reading background flux from file
    // Load and subtract off the background fields
    string loadFname = "dataPlane/"+fluxXFname+".h5";
    // Check if the file can be accessed by trying to open it
    ifstream in(loadFname.c_str(), ios::binary);
    if ( !in.good() )
    {
      cout << "File " << loadFname << " cannot be accessed. Run a background run first. bkg has to enter in the outdir path...\n";
      return 1;
    }
    in.close();

    field.set_output_directory("dataPlane");
    fluxYReflected.load_hdf5(field, fluxXFname.c_str());  
    fluxYReflected.scale_dfts(-1.0);
  }
  field.set_output_directory(OUTDIR); 

  while ( field.time() < tEnd )
  {
    field.step();

    // Get field amplitude
    complex<double> fieldAmpTrans = field.get_field(fieldComp, meep::vec(XSIZE/2.0, fluxPlanePosY));
    complex<double> fieldAmpRefl = field.get_field(fieldComp, meep::vec(XSIZE/2.0, fluxRefPlanePosY));
    fieldTransmitted.push_back(real(fieldAmpTrans));
    fieldReflection.push_back(real( fieldAmpRefl ));
    timepoints.push_back(field.time());

    #ifdef OUTPUT_HDF5
      if ( field.time() > nextOutputTime )
      {
        field.output_hdf5(fieldComp, vol.surroundings());
        nextOutputTime += dt;
      }
    #endif
  } 
  #ifdef OUTPUT_HDF5
    field.output_hdf5(fieldComp, vol.surroundings());
  #endif

  if ( outdir.find("bkg") != string::npos )
  {
    // This is the background run --> save the fields 
    field.set_output_directory("dataPlane");
    fluxYReflected.save_hdf5(field, fluxXFname.c_str());
  }

  // Write transmitted flux to file
  string ddir(OUTDIR);
  string fluxOut("transmittedFlux.csv");
  fluxOut = ddir+"/"+fluxOut; 
  ofstream os(fluxOut.c_str());
  if ( !os.good() )
  {
    cout << "Problem when opening file " << fluxOut << endl;
  }
  else
  {
    os << "# Flux from transmitted wave\n";
    os << "# Volume description\n";
    os << "# Flux in y-direction\n";
    os << "# Corner 1: (x,y) = (" << dftVol.get_min_corner().x() << ","<<dftVol.get_min_corner().y() << ")\n";
    os << "# Corner 2: (x,y) = (" << dftVol.get_max_corner().x() << "," << dftVol.get_max_corner().y() << ")\n"; 
    os << "# EPS_HIGH = " << EPS_HIGH << endl;
    os << "# Frequency, Angle, Flux Y, Flux Reflection\n";
    double dfreq = fwidth/static_cast<double>(nfreq-1);
    double currentFreq = freq-fwidth/2.0;
    double *transmittedFlux = transFluxY.flux();
    double *reflectedFlux = fluxYReflected.flux();
    unsigned int numberOfNonPropagating = 0;
    for ( unsigned int i=0;i<nfreq;i++ )
    {
      double sinCurrentAngle = freq*sin(angle*PI/180.0)/currentFreq;
      double currentAngle = 0.0;
      if ( abs(sinCurrentAngle) > 1.0 )
      {
        numberOfNonPropagating++;
      }
      else
      {
        currentAngle = asin(sinCurrentAngle)*180.0/PI;
        os << currentFreq << "," << currentAngle << "," << transmittedFlux[i]/transYWidth;
        os << "," << reflectedFlux[i]/transYWidth << "\n";
      }
      currentFreq += dfreq;
    }
    os.close();
    delete [] transmittedFlux;
    delete [] reflectedFlux;
    cout << numberOfNonPropagating << " of " << nfreq << " modes do not propagate\n";
  }
  
  // Write monitor to file
  string monitorOut("realField.csv");
  monitorOut = ddir + "/" + monitorOut;
  os.open(monitorOut.c_str());
  if ( !os.good() )
  {
    cerr << "Error when opening file " << monitorOut << endl;
    return 1;
  }
  os << "# Field monitored at positions\n";
  os << "# Time, Real field transmitted, Real field Reflected (+incident)\n";
  for ( unsigned int i=0;i<timepoints.size();i++ )
  {
    os << timepoints[i] << "," << fieldTransmitted[i] << "," << fieldReflection[i] << "\n";
  }
  os.close(); 
  return 0;
}

