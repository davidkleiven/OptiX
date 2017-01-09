#ifndef CYL_WG_CONF_MAP_H
#define CYL_WG_CONF_MAP_H
#include "curvedWaveGuide2D.hpp"

/** Simulate a curved waveguide in uv-plane after a confirmal map */
class CurvedWGConfMap: public CurvedWaveGuideFD
{
public:
  CurvedWGConfMap(): CurvedWaveGuideFD("CurvedWGConfMap"){};

  /** Return the spatial varying refractive index. The waveguide occupies the region [-width,0.0] */
  virtual void getXrayMatProp( double u, double v, double &delta, double &beta ) const override;

  /** True if the position is inside the waveguide */
  virtual bool isInsideGuide( double u, double v ) const override;
protected:
  CurvedWGConfMap( const char* name ): CurvedWaveGuideFD(name){};

  /** Returns the start position of the waveguide */
  virtual double waveGuideStartX( double v ) const override;

  /** Returns the end position of the waveguide */
  virtual double waveGuideEndX( double v ) const override;
};
#endif
