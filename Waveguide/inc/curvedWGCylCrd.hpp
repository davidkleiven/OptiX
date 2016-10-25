#ifndef CYLINDRICAL_WAVEGUIDE_CYLCRD_H
#define CYLINDRICAL_WAVEGUIDE_CYLCRD_H
#include "curvedWaveGuide2D.hpp"
#include <jsoncpp/json/writer.h>
#include <complex>

typedef std::complex<double > cdouble;

class CurvedWGCylCrd: public CurvedWaveGuideFD
{
public:
  CurvedWGCylCrd(): CurvedWaveGuideFD("CurvedWGCylCrd"){};
  cdouble transverseBC( double z, WaveGuideFDSimulation::Boundary_t bnd ) const override;
  virtual void fillInfo( Json::Value &obj ) const override;
protected:
  bool isInsideGuide( double r, double theta ) const override;
  bool waveguideEnded( double r, double theta ) const override;
};
#endif
