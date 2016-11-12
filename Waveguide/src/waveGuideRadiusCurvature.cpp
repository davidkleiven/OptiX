#include "waveGuideRadiusCurvature.hpp"
#include "cladding.hpp"
#include <iostream>

using namespace std;

double WaveGuideLargeCurvature::potential( double x ) const
{
  if (( x > 0.0 ) || ( x < -width ))
  {
    //double claddingPotential = cladding->getPotential();
    //return claddingPotential - (2.0*wavenumber*wavenumber - 2.0*claddingPotential)*x/outerRadius;
    double delta = cladding->getDelta();
    return 2.0*delta*wavenumber*wavenumber -(1.0-2.0*delta)*2.0*x*wavenumber*wavenumber/outerRadius;
  }
  return -2.0*wavenumber*wavenumber*x/outerRadius;
}

void WaveGuideLargeCurvature::load( ControlFile &ctl )
{
  WaveGuide1DSimulation::load(ctl);
  outerRadius = ctl.get()["RadiusOfCurvature"].asDouble();
}
