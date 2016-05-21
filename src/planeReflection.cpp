#include <iostream>
#include "meep.hpp"
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>
#include <fstream>
#include <string>
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
const double YSIZE = 10.0;
const double SOURCE_Y = 8.0;
const double PI = acos(-1.0);
const complex<double> IMAG_UNIT(0,1.0);
const unsigned int NSTEPS = 20;


double dielectric(const meep::vec &pos)
{ 
  if ( pos.y() < 5.0 )
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
  EPS_HIGH = *argv[2];
  ANGLE = *argv[3];

  double resolution = 20.0; // pixels per distance

  // Initialize computational cell
  meep::grid_volume vol = meep::vol2d(XSIZE, YSIZE, resolution);

  // Add line source to get plane wave
  meep::vec srcCorner1(0.0, SOURCE_Y);
  meep::vec srcCorner2(XSIZE, SOURCE_Y);
  meep::volume srcvol(srcCorner1, srcCorner2);

  // Initalize structure. Add PML in y-direction
  //meep::structure srct(vol, dielectric, meep::pml(1.0, meep::Y));
  meep::structure srct(vol, dielectric, meep::pml(1.0));

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
  meep::volume dftVol = meep::volume(meep::vec(0.0,3.0), meep::vec(0.0,XSIZE));
  meep::volume dftVolX = meep::volume(meep::vec(XSIZE*0.75, 0.0), meep::vec(XSIZE*0.75, 4.0));
  unsigned int nfreq = 20;
  meep::dft_flux transFluxY = field.add_dft_flux_plane(dftVol, freq+fwidth, freq-fwidth, nfreq);
  meep::dft_flux transFluxX = field.add_dft_flux_plane(dftVolX, freq+fwidth, freq-fwidth, nfreq);
  

  unsigned int nOut = 20; // Number of output files
  //double dt = ( field.last_source_time() + NSTEPS )/static_cast<double>(nOut); // Timestep between output hdf5
  double dt = nfreq/(nOut*fwidth);
  double nextOutputTime = 0.0;

  // Put a field monitor at the center of the geometry 
  meep::vec monitorPos(XSIZE/2.0, 7.0);

  vector<double> fieldAtCenterReal; // Container for the real field component
  vector<double> timepoints; // Container for the timepoints

  // Main loop.
  // TODO: Check if NSTEPS is correct. Maybe tune for frequency resulution in DFT
  while ( field.time() < nfreq/fwidth )
  {
    field.step();

    // Get field amplitude
    complex<double> fieldAmp = field.get_field(meep::Ez, monitorPos);
    fieldAtCenterReal.push_back(real(fieldAmp));
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
  os << "# Flux from transmitted wave\n";
  os << "# Volume description\n";
  os << "# Flux in y-direction\n";
  os << "# Corner 1: (x,y) = (" << dftVol.get_min_corner().x() << ","<<dftVol.get_min_corner().y() << ")\n";
  os << "# Corner 2: (x,y) = (" << dftVol.get_max_corner().x() << "," << dftVol.get_max_corner().y() << ")\n"; 
  os << "# Flux in x-direction:\n";
  os << "# Corner 1: (x,y) = (" << dftVolX.get_min_corner().x() << ","<<dftVolX.get_min_corner().y() << ")\n";
  os << "# Corner 2: (x,y) = (" << dftVolX.get_max_corner().x() << "," << dftVolX.get_max_corner().y() << ")\n"; 
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
  

  // Write monitor to file
  string monitorOut("ezMonitor.csv");
  monitorOut = ddir + "/" + monitorOut;
  os.open(monitorOut.c_str());
  if ( !os.good() )
  {
    cout << "Problem when opening file " << monitorOut << endl;
    return 1;
  }
  os << "# Field monitored at positions\n";
  os << "# Time, Ez.real\n";
  for ( unsigned int i=0;i<timepoints.size();i++ )
  {
    os << timepoints[i] << "," << fieldAtCenterReal[i] << "\n";
  }
  os.close(); 
  return 0;
}

