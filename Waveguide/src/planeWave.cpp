#include "planeWave.hpp"
#include <cmath>

const double PI = acos(-1.0);
cdouble PlaneWave::get( double x, double z ) const
{
  cdouble im(0.0,1.0);
  return amplitude*exp(im*wavenumber*x*sin(angleWithZAxis)-im*wavenumber*z0);
}

void PlaneWave::setAngleDeg( double angle )
{
  angleWithZAxis = angle*PI/180.0;
}

void PlaneWave::info( Json::Value &obj ) const
{
  ParaxialSource::info( obj );
  obj["angleWithZAxisDeg"] = angleWithZAxis*180.0/PI;
}
