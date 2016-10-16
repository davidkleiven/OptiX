#ifndef STRAIGHT_WG_2D_H
#define STRAIGHT_WG_2D_H
#include "curvedWaveGuide2D.hpp"
#include <complex>
#include <vector>

class ControlFile;

typedef std::complex<double> cdouble;

// A straight waveguide can be view as a curved one with radius of curvature R = infty
class StraightWG2D: public CurvedWaveGuideFD
{
public:
  StraightWG2D(): CurvedWaveGuideFD("StraightWaveGuide"){};
  void fillInfo( Json::Value &obj ) const override final;
protected:
  bool isInsideGuide( double x, double z ) const override final;
  double waveGuideStartX( double z ) const override final;
  double waveGuideEndX( double z ) const override final;
};
#endif
