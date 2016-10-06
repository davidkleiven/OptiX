#include "waveGuideRadiusCurvature.hpp"
#include "cladding.hpp"

double WaveGuideLargeCurvature::potential( double x ) const
{
  double claddingPotential = cladding->getPotential();
  return claddingPotential - (2.0*wavenumber*wavenumber - 2.0*claddingPotential)*x/outerRadius;
}
