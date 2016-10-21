#include "planeWave.hpp"
#include <cmath>

cdouble PlaneWave::get( double x, double z ) const
{
  cdouble im(0.0,1.0);
  return amplitude*exp(-im*wavenumber*z0);
}
