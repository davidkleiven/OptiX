#ifndef WAVE_GUIDE_HUGE_CURVATURE_H
#define WAVE_GUIDE_HUGE_CURVATURE_H
#include "waveGuide.hpp"

class WaveGuideLargeCurvature: public WaveGuide1DSimulation
{
public:
  WaveGuideLargeCurvature():WaveGuide1DSimulation("LargeRadiusOfCurvature"){};
  void setRadiusOfCurvature( double R ){ outerRadius = R; };
  double potential( double x ) const override;
};
#endif
