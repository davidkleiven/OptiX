#include "linearRampWG.hpp"
#include "cladding.hpp"

double LinearRampWG::profile( double x, double z ) const
{
  double criteria = 2.0*x*R+z*z;

  bool isInside = (criteria > 0.0) && (criteria < 2.0*width*R);

  bool isInUpperLayer = (criteria > 2*width*R) && (criteria < 2.0*(1+widthFraction)*width*R );
  bool isInLowerLayer = ( criteria < 0.0 ) && (criteria > -2.0*widthFraction*width*R);
  double val = x+0.5*z*z/R;
  if ( isInside )
  {
    return 0.0;
  }
  else if ( isInUpperLayer )
  {
    return (val-width)/( widthFraction*width );
  }
  else if ( isInLowerLayer )
  {
    return -val/(widthFraction*width);
  }
  return 1.0;
}

void LinearRampWG::getXrayMatProp( double x, double z, double &delta, double &beta ) const
{
  delta = cladding->getDelta()*profile(x,z);
  beta = cladding->getBeta()*profile(x,z);
}

void LinearRampWG::fillInfo( Json::Value &obj ) const
{
  CurvedWaveGuideFD::fillInfo(obj);
  obj["widthFraction"] = widthFraction;
}
