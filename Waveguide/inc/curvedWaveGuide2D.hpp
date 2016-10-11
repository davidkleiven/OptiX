#ifndef WAVE_GUIDE_CURVED_2D_H
#define WAVE_GUIDE_CURVED_2D_H
#include "waveGuideFDSimulation.hpp"
#include <complex>

typedef std::complex<double> cdouble;

class CurvedWaveGuideFD: public WaveGuideFDSimulation
{
public:
  CurvedWaveGuideFD(): WaveGuideFDSimulation("CurvedWaveGuide2D"){};
  cdouble getRefractiveIndex( double x, double z ) const override final;
  void setBoundaryConditions() override final;
protected:
  double R;
  double width;
  bool isInsideGuide( double x, double z ) const;
};
#endif
