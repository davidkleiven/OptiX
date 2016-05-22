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
#include "dataToFile.h"
#define ODIRLEN 60

/**
* This program takes 3 command line arguments
* Arg 1 - out directory. Directory where the datafiles will be stored
* Arg 2 - Relative dielectriv permittivity in scattering region. The wave starts in a region with epsilon=1.
* Arg 3 - angle of incidence in degrees
*
* Example:
*   ./planeReflection data/planewave 2.25 30
*/
using namespace std;

const double EPS_LOW = 1.0;
double EPS_HIGH = 1.0;

double ANGLE = 0.0;
const double XSIZE = 20.0;
const double YSIZE = 20.0;
const double PML_THICK = 4.0;
const double SOURCE_Y = YSIZE-PML_THICK - 1.0;
const double YC_PLANE = YSIZE/2.0;
const double PI = acos(-1.0);
const complex<double> IMAG_UNIT(0,1.0);
const unsigned int NSTEPS = 20;


double dielectric(const meep::vec &pos)
{
  if ( pos.y() < (YC_PLANE + (pos.x()-XSIZE/2.0)*tan(ANGLE*PI/180.0)) )
  {
    return EPS_HIGH;
  }
  return EPS_LOW;
}

complex<double> amplitude(const meep::vec &pos)
{
  return 1.0;
}

//----------------- MAIN FUNCTION ----------------------------//
int main(int argc, char **argv)
{
  cout << "This program simulates scattering of a plane wave onto a smooth surface\n";
  meep::initialize mpi(argc, argv);

  // Read command line arguments
  if ( argc != 4 )
  {
    cout << "Usage: ./planeReflection.out <outi directory> <epsInScattered> <incidend angle>\n";
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
  ss >> ANGLE;

  // Check that angle is within range
  const double maxAngle = atan( 2.0*(YC_PLANE-PML_THICK)/XSIZE )*180.0/PI;
  if ( ANGLE > maxAngle )
  {
    cout << "The incident angle is too large\n";
    cout << "Maximum angle is " << maxAngle << endl;
    return 1;
  }
  else if ( ANGLE < 0.0 )
  {
    cout << "Negative angle given. Has to be in range [0,MAX_ANGLE)\n";
    return 1;
  }

  // Verify that the size of the domain is big enough (for debugging only)
  const double minHeight = 2.0*PML_THICK + 4.0;
  assert ( YSIZE > minHeight );
  
  double resolution = 20.0; // pixels per distance

  // Initialize computational cell
  meep::grid_volume vol = meep::vol2d(XSIZE, YSIZE, resolution);

  // Add line source to get plane wave
  meep::vec srcCorner1(0.0, SOURCE_Y);
  meep::vec srcCorner2(XSIZE, SOURCE_Y);
  meep::volume srcvol(srcCorner1, srcCorner2);

  // Initalize structure. Add PML in y-direction
  //meep::structure srct(vol, dielectric, meep::pml(1.0, meep::Y));
  meep::structure srct(vol, dielectric, meep::pml(PML_THICK));

  srct.Courant = 0.1;
  meep::fields field(&srct);
  
  // Add periodic boundary conditions in x-direction
  //double blochK = 2.0*PI/(XSIZE-1.0);
  //field.use_bloch( meep::X, blochK ); 

  field.set_output_directory(OUTDIR); 

  // Write dielectric function to file
  field.output_hdf5(meep::Dielectric, vol.surroundings()); 

  // Set source type. Use Gaussian, had some problems with the continous
  double freq = 0.5;
  double fwidth = 0.2;
  meep::gaussian_src_time src(freq, fwidth);
  field.add_volume_source(meep::Ez, src, srcvol, amplitude);

  // Add DFT fluxplane
  meep::volume dftVol = meep::volume(meep::vec(0.0,PML_THICK), meep::vec(XSIZE,PML_THICK));
  meep::volume dftVolX = meep::volume(meep::vec(XSIZE-PML_THICK-1.0, PML_THICK), meep::vec(XSIZE-PML_THICK-1.0, YC_PLANE));
  unsigned int nfreq = 20;
  meep::dft_flux transFluxY = field.add_dft_flux_plane(dftVol, freq+fwidth, freq-fwidth, nfreq);
  meep::dft_flux transFluxX = field.add_dft_flux_plane(dftVolX, freq+fwidth, freq-fwidth, nfreq);
  

  unsigned int nOut = 20; // Number of output files
  //double dt = ( field.last_source_time() + NSTEPS )/static_cast<double>(nOut); // Timestep between output hdf5
  double dt = nfreq/(nOut*fwidth);
  double nextOutputTime = 0.0;

  // Put a field monitor at the center of the geometry 
  meep::vec monitorPos(XSIZE/2.0, 7.0);
  meep::vec monitorTransPlanePos(XSIZE/2.0, PML_THICK);

  vector<double> fieldAtCenterReal; // Container for the real field component
  vector<double> timepoints; // Container for the timepoints
  vector<double> fieldAtFluxPlane; // Container for the real field component at the center of the flux plane

  // Main loop.
  // TODO: Check if NSTEPS is correct. Maybe tune for frequency resulution in DFT
  while ( field.time() < nfreq/fwidth )
  {
    field.step();

    // Get field amplitude
    complex<double> fieldAmp = field.get_field(meep::Ez, monitorPos);
    complex<double> fieldAmpTrans = field.get_field(meep::Ez, monitorTransPlanePos);
    fieldAtCenterReal.push_back(real(fieldAmp));
    fieldAtFluxPlane.push_back(real( fieldAmpTrans ));
    timepoints.push_back(field.time());

    if ( field.time() > nextOutputTime )
    {
      field.output_hdf5(meep::Ez, vol.surroundings());
      nextOutputTime += dt;
    }
  } 
  field.output_hdf5(meep::Ez, vol.surroundings());

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
    os << "# Flux in x-direction:\n";
    os << "# Corner 1: (x,y) = (" << dftVolX.get_min_corner().x() << ","<<dftVolX.get_min_corner().y() << ")\n";
    os << "# Corner 2: (x,y) = (" << dftVolX.get_max_corner().x() << "," << dftVolX.get_max_corner().y() << ")\n"; 
    os << "# EPS_HIGH = " << EPS_HIGH << endl;
    os << "# Frequency, Flux Y, Flux X\n";
    double dfreq = 2.0*fwidth/static_cast<double>(nfreq);
    double currentFreq = freq-fwidth;
    double *transmittedFlux = transFluxY.flux();
    double *transmittedFluxX = transFluxX.flux();
    for ( unsigned int i=0;i<nfreq;i++ )
    {
      os << currentFreq << "," << transmittedFlux[i] << "," << transmittedFluxX[i] << "\n";
      currentFreq += dfreq;
    }
    os.close();
    delete [] transmittedFlux;
    delete [] transmittedFluxX;
  }
  

  // Write monitor to file
  try
  {
    string monitorOut("ezMonitorSource.csv");
    monitorOut = ddir + "/" + monitorOut;
    monitorToFile(monitorOut, timepoints, fieldAtCenterReal, monitorPos);
    monitorOut = ddir + "/" + "ezMonitorTrans.csv";
    monitorToFile(monitorOut, timepoints, fieldAtFluxPlane, monitorTransPlanePos);
  }
  catch ( runtime_error &exc )
  {
    cout << exc.what() << endl;
  }
  catch (...)
  {
    cout << "An unexpected exception occured...";
  }
  return 0;
}

