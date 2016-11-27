#ifndef WAVEGUIDE_STRAIGHT_H
#define WAVEGUIDE_STRAIGHT_H
#include "waveGuide.hpp"

/** Class for 1D straight waveguide */
class WaveGuideStraight: public WaveGuide1DSimulation
{
public:
  WaveGuideStraight():WaveGuide1DSimulation("Straight"){};
  ~WaveGuideStraight();

  /** Get potential induced by the refractive index profile */
  double potential( double x ) const override;
};
#endif
