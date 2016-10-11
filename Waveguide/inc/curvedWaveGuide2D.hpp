#ifndef WAVE_GUIDE_CURVED_2D_H
#define WAVE_GUIDE_CURVED_2D_H
#include "waveGuideFDSimulation.hpp"
#include <complex>

typedef std::complex<double> cdouble;

class CurvedWaveGuideFD: public WaveGuideFDSimulation
{
public:
  CurvedWaveGuideFD(): WaveGuideFDSimulation("CurvedWaveGuide2D"){};
  void getXrayMatProp( double x, double z, double &delta, double &beta ) const override final;
  void setBoundaryConditions() override final;
  void fillInfo( Json::Value &obj ) const override final;
  void setRadiusOfCurvature( double newR ) { R = newR; };
  void setWidth( double newWidth ) { width = newWidth; };
protected:
  double R;
  double width;
  bool isInsideGuide( double x, double z ) const override final;
};
#endif
