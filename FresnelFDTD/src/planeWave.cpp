#include "planeWave.h"
#include <complex>
#include <cmath>

double PlaneWave::EPS_HIGH = 1.0;
double PlaneWave::KX = 1.0;

double PlaneWave::dielectric(const meep::vec &pos)
{
  if ( pos.y() < PlaneWave::YC_PLANE )
  {
    return PlaneWave::EPS_HIGH;
  }
  return PlaneWave::EPS_LOW;
}

std::complex<double> PlaneWave::amplitude(const meep::vec &pos)
{
  return exp(PlaneWave::IMAG_UNIT*PlaneWave::KX*pos.x());
}
