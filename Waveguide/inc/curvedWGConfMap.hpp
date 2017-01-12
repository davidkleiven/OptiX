#ifndef CYL_WG_CONF_MAP_H
#define CYL_WG_CONF_MAP_H
#include "curvedWaveGuide2D.hpp"

/** Simulate a curved waveguide in uv-plane after a confirmal map */
class CurvedWGConfMap: public CurvedWaveGuideFD
{
public:
  enum class Curvature_t{ CONCAVE, CONVEX };
  CurvedWGConfMap(): CurvedWaveGuideFD("CurvedWGConfMap"){};

  /** Get curvature */
  Curvature_t getCurvature() const { return curvature; };

  /** Set curvature of the waveguide */
  void setCurvature( Curvature_t curv ){ curvature = curv; };

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

  Curvature_t curvature{Curvature_t::CONCAVE};
};
#endif
