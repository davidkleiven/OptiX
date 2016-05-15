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
const double EPS_HIGH = 1.5;
char OUTDIR[ODIRLEN] = "dataPlane/NormalInc";
const double XSIZE = 10.0;
const double YSIZE = 10.0;
const double SOURCE_Y = 7.0;
const double ANGLE=10.0;
const double PI = acos(-1.0);
const complex<double> IMAG_UNIT(0,1.0);
const string OUT_MONITOR_FNAME("dataPlane/NormalInc/ezMonitor.csv");


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

double sourceX()
{
  double xcenter = XSIZE/2.0;
  double delta = SOURCE_Y*tan(ANGLE*PI/180.0);
  if ( delta > xcenter )
  {
    cout << "The angle is too large. Increase size in xdirection\n";
    return xcenter;
  }
  return xcenter-delta;
} 

int main(int argc, char **argv)
{
  cout << "This program simulates scattering of a plane wave onto a smooth surface\n";
  meep::initialize mpi(argc, argv);

  double resolution = 20.0; // pixels per distance
  meep::grid_volume vol = meep::vol2d(XSIZE, YSIZE, resolution);
  meep::vec srcCorner1(0.0, SOURCE_Y);
  meep::vec srcCorner2(XSIZE, SOURCE_Y);
  meep::volume srcvol(srcCorner1, srcCorner2);

  meep::structure srct(vol, dielectric, meep::pml(1.0, meep::Y));
  srct.Courant = 0.4;
  meep::fields field(&srct);
  double blochK = 2.0*PI/(XSIZE-1.0);
  field.use_bloch( meep::X, blochK ); 

  char* outdirheap = new char[ODIRLEN]; // This is deleted in the destructor of the meep::fields class
  memcpy(outdirheap, OUTDIR, ODIRLEN);
  
  // Check that field.outdir is not allocated
  if ( field.outdir != NULL )
  {
    cout << "The out directory is already allocated! Will deallocate it!\n";
    cout << "Content:\n";
    unsigned int indx=0;
    while ( (field.outdir[indx] != '\0') && (indx < 100000) )
    {
      cout << field.outdir[indx];
      indx++;
    }
    cout << endl;
    delete [] field.outdir;
  }
  field.outdir = outdirheap;
  field.output_hdf5(meep::Dielectric, vol.surroundings()); 

  double freq = 0.3;
  double freqwidth = 0.1;
  meep::gaussian_src_time src(freq, freqwidth);
  field.add_volume_source(meep::Ez, src, srcvol, amplitude);

  unsigned int nOut = 20;
  double dt = field.last_source_time()/static_cast<double>(nOut);
  double nextOutputTime = 0.0;
  meep::vec monitorPos(XSIZE/2.0, 6.0);
  vector<double> fieldAtCenterReal;
  vector<double> fieldAtCenterImag;
  vector<double> timepoints;
  while ( field.time() < field.last_source_time() )
  {
    field.step();
    complex<double> fieldAmp = field.get_field(meep::Ez, monitorPos);
    fieldAtCenterReal.push_back(real(fieldAmp));
    fieldAtCenterImag.push_back(imag(fieldAmp));
    timepoints.push_back(field.time());
    if ( field.time() > nextOutputTime )
    {
      field.output_hdf5(meep::Ez, vol.surroundings());
      nextOutputTime += dt;
    }
  } 
  field.output_hdf5(meep::Ez, vol.surroundings());


  // Write monitor to file
  ofstream os(OUT_MONITOR_FNAME.c_str());
  if ( !os.good() )
  {
    cout << "Problem when opening file " << OUT_MONITOR_FNAME << endl;
    return 1;
  }

  os << "# Field monitored at positions\n";
  os << "# Time, Ez.real, Ez.imag\n";
  for ( unsigned int i=0;i<timepoints.size();i++ )
  {
    os << timepoints[i] << "," << fieldAtCenterReal[i] << "," << fieldAtCenterImag[i] << "\n";
  }
  os.close(); 
  return 0;
}

