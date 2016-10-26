#ifndef GAUSSIAN_PROFILE_WG
#define GAUSSIAN_PROFILE_WG
#include "curvedWaveGuide2D.hpp"

class GaussianWG: public CurvedWaveGuideFD
{
public:
  GaussianWG(): CurvedWaveGuideFD("curvedWGGaussianProfile"){};
  virtual void getXrayMatProp( double x, double z, double &delta, double &beta ) const override;
private:
  double profile( double x, double z ) const;
};

#endif
