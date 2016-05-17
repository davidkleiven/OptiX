#include <iostream>
#include "meep.hpp"
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>
#include <fstream>
#include <string>
#define ODIRLEN 60

using namespace std;

const double EPS_LOW = 1.0;
const double EPS_HIGH = 2.25;
char OUTDIR[ODIRLEN] = "dataPlane/NormalInc";

const double XSIZE = 10.0;
const double YSIZE = 10.0;
const double SOURCE_Y = 8.0;
const double ANGLE=10.0;
const double PI = acos(-1.0);
const complex<double> IMAG_UNIT(0,1.0);
const unsigned int NSTEPS = 20;


double dielectric(const meep::vec &pos)
{ 
  if ( pos.y() < 3.0 )
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

  double resolution = 10.0; // pixels per distance

  // Initialize computational cell
  meep::grid_volume vol = meep::vol2d(XSIZE, YSIZE, resolution);

  // Add line source to get plane wave
  meep::vec srcCorner1(0.0, SOURCE_Y);
  meep::vec srcCorner2(XSIZE, SOURCE_Y);
  meep::volume srcvol(srcCorner1, srcCorner2);

  // Initalize structure. Add PML in y-direction
  meep::structure srct(vol, dielectric, meep::pml(1.0, meep::Y));

  srct.Courant = 0.1;
  meep::fields field(&srct);
  
  // Add periodic boundary conditions in x-direction
  double blochK = 2.0*PI/(XSIZE-1.0);
  field.use_bloch( meep::X, blochK ); 


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
  unsigned int nfreq = 20;
  meep::dft_flux transFlux = field.add_dft_flux_plane(dftVol, freq+fwidth, freq-fwidth, nfreq);

  unsigned int nOut = 20; // Number of output files
  double dt = ( field.last_source_time() + NSTEPS )/static_cast<double>(nOut); // Timestep between output hdf5
  double nextOutputTime = 0.0;

  // Put a field monitor at the center of the geometry 
  meep::vec monitorPos(XSIZE/2.0, 7.0);

  vector<double> fieldAtCenterReal; // Container for the real field component
  vector<double> timepoints; // Container for the timepoints

  // Main loop.
  // TODO: Check if NSTEPS is correct. Maybe tune for frequency resulution in DFT
  while ( field.time() < nfreq/fwidth + NSTEPS )
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
  os << "# Corner 1: (x,y) = (" << dftVol.get_min_corner().x() << ","<<dftVol.get_min_corner().y() << ")\n";
  os << "# Corner 2: (x,y) = (" << dftVol.get_max_corner().x() << "," << dftVol.get_max_corner().y() << ")\n"; 
  double dfreq = 2.0*fwidth/static_cast<double>(nfreq);
  double currentFreq = freq-fwidth;
  double *transmittedFlux = transFlux.flux();
  for ( unsigned int i=0;i<nfreq;i++ )
  {
    os << currentFreq << "," << transmittedFlux[i] << "\n";
    currentFreq += dfreq;
  }
  os.close();
  

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

