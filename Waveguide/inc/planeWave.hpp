#ifndef PLANE_WAVE_SOURCE_H
#define PLANE_WAVE_SOURCE_H
#include "paraxialSource.hpp"
#include <jsoncpp/json/writer.h>

class PlaneWave: public ParaxialSource
{
public:
  PlaneWave(): ParaxialSource("planeWave"){};
  cdouble get( double x, double z ) const override final;
  void setAngleDeg( double angle );
  void info( Json::Value &obj ) const override;
private:
  double angleWithZAxis{0.0};
};
#endif
