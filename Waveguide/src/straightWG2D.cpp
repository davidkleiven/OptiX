#include "straightWG2D.hpp"

bool StraightWG2D::isInsideGuide( double x, double z ) const
{
  return ( x > 0.0 ) && ( x < width );
}

void StraightWG2D::fillInfo( Json::Value &obj ) const
{
    obj["Width"] = width;
}

double StraightWG2D::waveGuideStartX( double z ) const
{
  return 0.0;
}

double StraightWG2D::waveGuideEndX( double z ) const
{
  return width;
}
