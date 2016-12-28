#ifndef PLANE_WAVE_SOURCE_H
#define PLANE_WAVE_SOURCE_H
#include "paraxialSource.hpp"
#include <json/writer.h>

/** Class for handling the plane wave source */
class PlaneWave: public ParaxialSource
{
public:
  PlaneWave(): ParaxialSource("planeWave"){};

  /** Evaluate the field at position x, z */
  cdouble get( double x, double z ) const override final;

  /** Set angle in degrees with the z-axis */
  void setAngleDeg( double angle );

  /** Fill JSON object with parameters specifici to this class */
  void info( Json::Value &obj ) const override;
private:
  double angleWithZAxis{0.0};
};
#endif
