#ifndef WAVE_GUIDE_CURVED_2D_H
#define WAVE_GUIDE_CURVED_2D_H
#include "waveGuideFDSimulation.hpp"
#include <complex>
#include <vector>

class ControlFile;

typedef std::complex<double> cdouble;

class CurvedWaveGuideFD: public WaveGuideFDSimulation
{
public:
  CurvedWaveGuideFD(): WaveGuideFDSimulation("CurvedWaveGuide2D"){};
  void getXrayMatProp( double x, double z, double &delta, double &beta ) const override final;
  void setBoundaryConditions() override final;
  void setRadiusOfCurvature( double newR ) { R = newR; };
  void setWidth( double newWidth ) { width = newWidth; };
  cdouble transverseBC( double z ) const override final;
  void computeTransmission( double step );
  void saveTransmission( ControlFile &ctl ) const;

  // Virtual functions
  virtual void fillInfo( Json::Value &obj ) const override;
  virtual void init( const ControlFile &ctl ) override;
protected:
  CurvedWaveGuideFD( const char *name): WaveGuideFDSimulation(name){};
  double R;
  double width;
  virtual bool isInsideGuide( double x, double z ) const override;
  virtual double waveGuideStartX( double z ) const;
  virtual double waveGuideEndX( double z ) const;
  std::vector<double> transmission;
  double stepWhenComputingTransmission{0.0};
};
#endif
