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
