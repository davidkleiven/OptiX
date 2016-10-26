#ifndef LINEAR_RAMP_WG_H
#define LINEAR_RAMP_WG_H
#include "curvedWaveGuide2D.hpp"

class LinearRampWG: public CurvedWaveGuideFD
{
public:
  LinearRampWG(): CurvedWaveGuideFD("CurvedWGLinearRamp"){};
  void setWidthFraction( double f ){ widthFraction=f; };
  virtual void getXrayMatProp( double x, double z, double &delta, double &beta ) const override;
  virtual void fillInfo( Json::Value &obj ) const override;
private:
  double widthFraction{0.0};
  double profile( double x, double z ) const;
};
#endif
