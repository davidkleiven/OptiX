#include "cylindricalParaxialEquation.hpp"

double CylindricalParaxialEquation::F( double r, double theta ) const
{
  return (1.0/R);
}

double CylindricalParaxialEquation::G( double r, double theta ) const
{
  //return 1.0;
  return (1.0/R)*(1.0-r/R);
  return 1.0/r;
}

double CylindricalParaxialEquation::H( double r, double theta ) const
{
  return (R+r);
  return 1.0;
  return r;
}

double CylindricalParaxialEquation::J( double r, double theta ) const
{
  return r/R;
}
