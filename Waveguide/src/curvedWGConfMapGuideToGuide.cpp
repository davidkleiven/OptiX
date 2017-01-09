#include "curvedWGConfMapGuideToGuide.hpp"

void CurvedWGConfMapGuideToGuide::getXrayMatProp( double u, double v, double &delta, double &beta ) const
{
  double ratio = R/rTarget;
  double transformedU = sign*u*R/rTarget - sign*R*log(R/rTarget);
  CurvedWGConfMap::getXrayMatProp( transformedU, v, delta, beta );
  delta *= pow(ratio,2);
  beta *= pow(ratio,2);
}

void CurvedWGConfMapGuideToGuide::setSign( int newsign )
{
  sign = newsign > 0 ? 1:-1;
}
