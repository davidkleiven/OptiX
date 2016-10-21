#include "gaussianBeam.hpp"
#include <cmath>

using namespace std;
const double PI = acos(-1.0);

double GaussianBeam::rayleighRange() const
{
  return PI*waist*waist/getWavelength();
}

double GaussianBeam::guoyPhase( double z ) const
{
  return atan( z/rayleighRange() );
}

double GaussianBeam::inverseRadiusOfCurvature( double z ) const
{
  double zR = rayleighRange();
  return z/( z*z + zR*zR);
}

double GaussianBeam::spotSize( double z ) const
{
  double zR = rayleighRange();
  return waist*sqrt( 1.0 + pow(z/zR,2));
}

cdouble GaussianBeam::get( double x, double z ) const
{
  cdouble im(0.0,1.0);
  double waistRatio = waist/spotSize( z-z0 );
  double gaussianFactor = exp(-pow(x/spotSize(z-z0),2) );
  cdouble phaseFactor = exp(-im*wavenumber*z0 + 0.5*im*wavenumber*x*x*inverseRadiusOfCurvature(z-z0) \
                          -im*guoyPhase(z-z0));
  return amplitude*waistRatio*gaussianFactor*phaseFactor;
}

double GaussianBeam::beamDivergence() const
{
  return getWavelength()/(PI*waist);
}
