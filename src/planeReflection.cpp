#include <iostream>
#include "meep.hpp"

using namespace std;

const double EPS_LOW = 1.0;
const double EPS_HIGH = 1.5;

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

  field.output_hdf5(meep::Dielectric, vol.surroundings()); 
  return 0;
}

