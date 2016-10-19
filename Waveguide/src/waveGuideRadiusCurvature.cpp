#include "waveGuideRadiusCurvature.hpp"
#include "cladding.hpp"
#include <iostream>

using namespace std;

double WaveGuideLargeCurvature::potential( double x ) const
{
  if (( x > 0.0 ) || ( x < -width ))
  {
    double claddingPotential = cladding->getPotential();
    return claddingPotential - (2.0*wavenumber*wavenumber - 2.0*claddingPotential)*x/outerRadius;
  }
  return -2.0*wavenumber*wavenumber*x/outerRadius;
}

void WaveGuideLargeCurvature::load( ControlFile &ctl )
{
  WaveGuide1DSimulation::load(ctl);
  outerRadius = ctl.get()["RadiusOfCurvature"].asDouble();
}
