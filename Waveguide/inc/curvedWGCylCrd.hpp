#ifndef CYLINDRICAL_WAVEGUIDE_CYLCRD_H
#define CYLINDRICAL_WAVEGUIDE_CYLCRD_H
#include "curvedWaveGuide2D.hpp"
#include <json/writer.h>
#include <complex>

typedef std::complex<double > cdouble;

/** Class for handling the curved waveguide in cylindrical coordinates */
class CurvedWGCylCrd: public CurvedWaveGuideFD
{
public:
  CurvedWGCylCrd(): CurvedWaveGuideFD("CurvedWGCylCrd"){};

  /** Returns the boundary condition at r=rmin and r=rmax (assumed to be the same on both locations)*/
  cdouble transverseBC( double z, WaveGuideFDSimulation::Boundary_t bnd ) const override;

  /** Fill JSON object with information relevant for this geometry */
  virtual void fillInfo( Json::Value &obj ) const override;
protected:
  bool isInsideGuide( double r, double theta ) const override;
  bool waveguideEnded( double r, double theta ) const override;

  /** Returns the lower waveguide border at positon z. Required when computing the transmission */
  double waveGuideStartX( double z ) const override;

  /** Reuturns the upper waveguide border at position z. Required when computing the transmission */
  double waveGuideEndX( double z ) const override;
};
#endif
