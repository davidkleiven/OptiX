#ifndef CURVED_WG_CONF_MAP_POSITIVE_H
#define CURVED_WG_CONF_MAP_POSITIVE_H
#include "curvedWGConfMap.hpp"

/** Conformal map from a waveguide in uv-plane with a certain radius to another with a different radius */
class CurvedWGConfMapGuideToGuide: public CurvedWGConfMap
{
public:
  CurvedWGConfMapGuideToGuide(): CurvedWGConfMap("CurvedWGConfMapGuideToGuide"){};

  /** Set the radius of the target coordinate system */
  void setTargetRadius( double r ){ rTarget = r; };

  /** Set the sign of the curvature */
  void setSign( int sign );

  /** Return the spatial varying refractive index. The waveguide occupies the region [-width,0.0] */
  virtual void getXrayMatProp( double u, double v, double &delta, double &beta ) const override;
private:
  int sign{1};
  double rTarget{40.0E6};
};

#endif
