#ifndef WG_COMPLEX_POTENTIAL_H
#define WG_COMPLEX_POTENTIAL_H
#include "waveGuideRadiusOfCuvature.hpp"

/** Class handling the 1D schrodinger equation witha complex potential */
class WaveguideComplex: public WaveGuideLargeCurvature
{
  WaveguideComplex();

  /** Get the complex values potential */
  cdouble complexPotential( double x ) const override;

  /** Run the simulation */
  cdouble solve() override;
};

#endif
