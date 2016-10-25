#include "beamSplitter.hpp"
#include <cmath>

using namespace std;

void BeamSplitter::fillInfo( Json::Value &obj ) const
{
  obj["anlgeDeg"] = angleWithZAxisDeg;
  obj["splitStart"] = splitStart;
  obj["width"] = width;
}

bool BeamSplitter::isInsideGuide( double x, double z ) const
{
  const double PI = acos(-1.0);
  if ( z > splitStart )
  {
    double slope = angleWithZAxisDeg*PI/180.0;
    double xmin = slope*(z-splitStart);
    double xmax = xmin+width;
    return (abs(x) > xmin) && (abs(x) < xmax);
  }
  return abs(x) < width;
}
