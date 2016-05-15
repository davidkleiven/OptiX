#include <iostream>
#include "meep.hpp"
#include <cstring>
#include <cmath>
#include <complex>

using namespace std;

const double EPS_LOW = 1.0;
const double EPS_HIGH = 1.5;
char OUTDIR[5] = "data";
const double XSIZE = 5.0;
const double YSIZE = 10.0;
const double SOURCE_Y = 9.0;
const double ANGLE=10.0;
const double PI = acos(-1.0);
const complex<double> IMAG_UNIT(0,1.0);


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
  //return 1.0;
  double lambda = 1.0/0.3;
  double kx = 2.0*PI/lambda;
  return exp(IMAG_UNIT*kx*pos.x());
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
  meep::vec srcCorner1(0.0, 9.0);
  meep::vec srcCorner2(XSIZE, 9.0);
  meep::volume srcvol(srcCorner1, srcCorner2);

  meep::structure srct(vol, dielectric, meep::pml(1.0));
  meep::fields field(&srct);

  char* outdirheap = new char[5]; // This is deleted in the destructor of the meep::fields class
  memcpy(outdirheap, OUTDIR, 5);
  
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
  //meep::vec srcpos(sourceX(), SOURCE_Y);
  //field.add_point_source(meep::Ez, src, srcpos);
  field.add_volume_source(meep::Ez, src, srcvol, amplitude);

  unsigned int nOut = 40;
  double dt = field.last_source_time()/static_cast<double>(nOut);
  double nextOutputTime = 0.0;
  while ( field.time() < field.last_source_time() )
  {
    field.step();
    if ( field.time() > nextOutputTime )
    {
      field.output_hdf5(meep::Ez, vol.surroundings());
      nextOutputTime += dt;
    }
  } 
  field.output_hdf5(meep::Ez, vol.surroundings());
  return 0;
}

