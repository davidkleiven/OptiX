#include "curvedWGCylCrd.hpp"
#include <cmath>
#include <iostream>

using namespace std;

bool CurvedWGCylCrd::isInsideGuide( double r, double theta ) const
{
  return ((r > 0.0) && (r < width));
}

bool CurvedWGCylCrd::waveguideEnded( double r, double theta ) const
{
  return WaveGuideFDSimulation::waveguideEnded( r, R*theta );
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

cdouble CurvedWGCylCrd::padExitField( double r, double theta ) const
{
  return WaveGuideFDSimulation::padExitField( r, R*theta );
}
