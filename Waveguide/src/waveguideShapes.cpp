#include "waveguideShapes.hpp"

bool SquareWell::isInside( double x, double y, double z ) const
{
  bool yInside = ( y > 0.0 );
  bool xInside = (( x > -width/2.0 ) && ( x < width/2.0 ));
  bool yAboveTop = ( y > depth );

  return (yInside && xInside) || yAboveTop;
}
