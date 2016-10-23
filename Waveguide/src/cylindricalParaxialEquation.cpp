#include "cylindricalParaxialEquation.hpp"

double CylindricalParaxialEquation::F( double r, double theta ) const
{
  return 1.0/r;
}

double CylindricalParaxialEquation::G( double r, double theta ) const
{
  return 1.0/r;
}

double CylindricalParaxialEquation::H( double r, double theta ) const
{
  return r;
}
