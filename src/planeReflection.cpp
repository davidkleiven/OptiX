#include <iostream>
#include "meep.hpp"
#include <cstring>

using namespace std;

const double EPS_LOW = 1.0;
const double EPS_HIGH = 1.5;
char OUTDIR[5] = "data";


double dielectric(const meep::vec &pos)
{ 
  if ( pos.y() > 5.0 )
  {
    return EPS_HIGH;
  }
  return EPS_LOW;
}

int main(int argc, char **argv)
{
  cout << "This program simulates scattering of a plane wave onto a smooth surface\n";
  meep::initialize mpi(argc, argv);

  double resolution = 20.0; // pixels per distance
  double xsize = 5.0;
  double ysize = 10.0;
  meep::grid_volume vol = meep::vol2d(xsize, ysize, resolution);
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
  meep::vec srcpos(1.0, 8.0);
  field.add_point_source(meep::Ez, src, srcpos);

  unsigned int nOut = 4;
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

