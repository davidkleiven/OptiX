#include "paraxialEquation.hpp"

double ParaxialEquation::F( double x, double z ) const
{
  return 1.0;
}

double ParaxialEquation::G( double x, double z ) const
{
  return 1.0;
}

double ParaxialEquation::H( double x, double z ) const
{
  return 1.0;
}

double ParaxialEquation::J( double x, double z ) const
{
  return 0.0;
}

cdouble ParaxialEquation::phaseFactor( double k, double z ) const
{
  cdouble im(0.0,1.0);
  return exp(im*k*z);
}
