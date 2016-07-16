#ifndef PLANE_WAVE_H
#define PLANE_WAVE_H
#include "meep.hpp"
#include <complex>

namespace PlaneWave{
  const double EPS_LOW = 1.0;
  extern double EPS_HIGH;

  const double XSIZE = 5.0;
  const double YSIZE = 18.0;
  const double PML_THICK = 3.0;
  const double SOURCE_Y = YSIZE - PML_THICK - 1.0;
  const double YC_PLANE = YSIZE/2.0;
  //const double YC_PLANE = PML_THICK+1.0;
  const double PI = acos(-1.0);
  const std::complex<double> IMAG_UNIT(0,1.0);
  const unsigned int NSTEPS = 20;
  extern double KX;

  double dielectric(const meep::vec &pos);

  std::complex<double> amplitude(const meep::vec &pos);
};

#endif
