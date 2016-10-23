#ifndef PARAXIAL_EQUATION_H
#define PARAXIAL_EQUATION_H
#include <complex>

typedef std::complex<double> cdouble;

/**
* Interface class for the following equation
*
* F(x,z)dA/dz = a*G(x,z)d/dx (H(x,z)dA/dx) - K*A + i*k*J(x,z)*A
*
* a and K are related to the material parameters and are handles elsewhere
* k is the wavenumber
* Default is cartesian coordinates where F = G = H = 1
*/
class ParaxialEquation
{
public:
  ParaxialEquation(){};
  virtual double F( double x, double z ) const;
  virtual double G( double x, double z ) const;
  virtual double H( double x, double z ) const;
  virtual double J( double x, double z ) const;
};
#endif
