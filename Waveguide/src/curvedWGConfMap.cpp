#include "curvedWGConfMap.hpp"
#include "cladding.hpp"
#include <iostream>

using namespace std;

void CurvedWGConfMap::getXrayMatProp( double u, double v, double &delta, double &beta ) const
{
  delta = -u/R;
  beta = 0.0;
  if (!isInsideGuide(u,v))
  {
    delta += cladding->getDelta();
    beta = cladding->getBeta();
  }
  return;
}

double CurvedWGConfMap::waveGuideStartX( double v ) const
{
  return -width;
}

double CurvedWGConfMap::waveGuideEndX( double v ) const
{
  return 0.0;
}

bool CurvedWGConfMap::isInsideGuide( double u, double v ) const
{
  return ( u > -width ) && ( u < 0.0 );
}
