#include "cylindricalParaxialEquation.hpp"

cdouble CylindricalParaxialEquation::F( double r, double theta ) const
{
  return 1.0/r;
}

cdouble CylindricalParaxialEquation::G( double r, double theta ) const
{
  return 1.0/r;
}

cdouble CylindricalParaxialEquation::H( double r, double theta ) const
{
  return r;
}
