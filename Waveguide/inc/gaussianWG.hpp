#ifndef GAUSSIAN_PROFILE_WG
#define GAUSSIAN_PROFILE_WG
#include "curvedWaveGuide2D.hpp"

/** Gaussian refractive index profile */
class GaussianWG: public CurvedWaveGuideFD
{
public:
  GaussianWG(): CurvedWaveGuideFD("curvedWGGaussianProfile"){};

  /** Return the refractive index according to a Gaussian weighting factor*/
  virtual void getXrayMatProp( double x, double z, double &delta, double &beta ) const override;
private:
  /** Return the Gaussian weighting factor */
  double profile( double x, double z ) const;
};

#endif
