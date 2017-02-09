#include "waveguide3DPaths.hpp"

void StraightPath::get( double z, double &x, double &y ) const
{
  x = 0.0;
  y = 0.0;
}

void ParabolicPath::get( double z, double &x, double &y ) const
{
  y = 0.0;
  x = -z*z/(2.0*R*1E6);
}
