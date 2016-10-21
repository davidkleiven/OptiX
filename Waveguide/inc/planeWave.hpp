#ifndef PLANE_WAVE_SOURCE_H
#define PLANE_WAVE_SOURCE_H
#include "paraxialSource.hpp"

class PlaneWave: public ParaxialSource
{
public:
  PlaneWave(): ParaxialSource("planeWave"){};
  cdouble get( double x, double z ) const override final;
};
#endif
