#ifndef WAVEGUIDE_STRAIGHT_H
#define WAVEGUIDE_STRAIGHT_H
#include "waveGuide.hpp"

class WaveGuideStraight: public WaveGuide1DSimulation
{
public:
  WaveGuideStraight():WaveGuide1DSimulation("Straight"){};
  ~WaveGuideStraight();
  double potential( double x ) const override;
};
#endif
