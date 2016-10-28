#include "curvedWGCylCrd.hpp"
#include <cmath>

bool CurvedWGCylCrd::isInsideGuide( double r, double theta ) const
{
  return ((r > 0) && (r < width));
}

bool CurvedWGCylCrd::waveguideEnded( double r, double theta ) const
{
  return WaveGuideFDSimulation::waveguideEnded( r, R*theta );
}

cdouble CurvedWGCylCrd::transverseBC( double theta, WaveGuideFDSimulation::Boundary_t bnd ) const
{
  double r = 0.0;
  if ( bnd == WaveGuideFDSimulation::Boundary_t::TOP )
  {
    r = R+xDisc->max;
  }
  else
  {
    r = R+xDisc->min;
  }
  return WaveGuideFDSimulation::transverseBC( r*sin(theta), bnd );
}

void CurvedWGCylCrd::fillInfo( Json::Value &obj ) const
{
  CurvedWaveGuideFD::fillInfo(obj);
  obj["crd"]="cylindrical";
}

double CurvedWGCylCrd::waveGuideStartX( double z ) const
{
  return 0.0;
}

double CurvedWGCylCrd::waveGuideEndX( double z ) const
{
  return width;
}
