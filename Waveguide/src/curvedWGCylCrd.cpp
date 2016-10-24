#include "curvedWGCylCrd.hpp"

bool CurvedWGCylCrd::isInsideGuide( double r, double theta ) const
{
  return (r > 0) && (r < width);
}

cdouble CurvedWGCylCrd::transverseBC( double theta, WaveGuideFDSimulation::Boundary_t bnd ) const
{
  return WaveGuideFDSimulation::transverseBC( R*theta, bnd );
}

void CurvedWGCylCrd::fillInfo( Json::Value &obj ) const
{
  CurvedWaveGuideFD::fillInfo(obj);
  obj["crd"]="cylindrical";
}
