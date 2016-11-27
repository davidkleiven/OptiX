#ifndef LINEAR_RAMP_WG_H
#define LINEAR_RAMP_WG_H
#include "curvedWaveGuide2D.hpp"

/** Class for handling a linear refractive index profile */
class LinearRampWG: public CurvedWaveGuideFD
{
public:
  LinearRampWG(): CurvedWaveGuideFD("CurvedWGLinearRamp"){};

  /** Set the width over which the linear profile changes from 0 to 1 given in units of the waveguide width */
  void setWidthFraction( double f ){ widthFraction=f; };

  /** Get the refractive index according to the linear weighting factor */
  virtual void getXrayMatProp( double x, double z, double &delta, double &beta ) const override;

  /** Fill JSON object with parameters specific to this class */
  virtual void fillInfo( Json::Value &obj ) const override;
private:
  double widthFraction{0.0};

  /** Return the linear profil weighting factor */
  double profile( double x, double z ) const;
};
#endif
