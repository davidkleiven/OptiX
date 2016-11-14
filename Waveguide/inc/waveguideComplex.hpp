#ifndef WG_COMPLEX_POTENTIAL_H
#define WG_COMPLEX_POTENTIAL_H
#include "waveGuideRadiusOfCuvature.hpp"

class WaveguideComplex: public WaveGuideLargeCurvature
{
  WaveguideComplex();
  cdouble complexPotential( double x ) const override;
  cdouble solve() override;
};

#endif
