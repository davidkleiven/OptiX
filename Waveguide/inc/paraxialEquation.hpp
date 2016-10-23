#ifndef PARAXIAL_EQUATION_H
#define PARAXIAL_EQUATION_H
#include <complex>

typedef std::complex<double> cdouble;

/**
* Interface class for the following equation
*
* F(x,z)dA/dz = a*G(x,z)d/dx (H(x,z)dA/dx) - K*A
* Default is cartesian coordinates where F = G = H = 1
*/
class ParaxialEquation
{
public:
  ParaxialEquation(){};
  virtual cdouble F( double x, double z ) const;
  virtual cdouble G( double x, double z ) const;
  virtual cdouble H( double x, double z ) const;
};
#endif
