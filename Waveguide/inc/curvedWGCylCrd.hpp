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

  /** Fill JSON object with information relevant for this geometry */
  virtual void fillInfo( Json::Value &obj ) const override;
protected:
  bool isInsideGuide( double r, double theta ) const override;
  bool waveguideEnded( double r, double theta ) const override;

  /** Returns the lower waveguide border at positon z. Required when computing the transmission */
  double waveGuideStartX( double z ) const override;

  /** Reuturns the upper waveguide border at position z. Required when computing the transmission */
  double waveGuideEndX( double z ) const override;

  /** Pad the signal corresponding to the unaffected source */
  virtual cdouble padExitField( double r, double theta ) const override;
};
#endif
