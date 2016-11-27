#ifndef STRAIGHT_WG_2D_H
#define STRAIGHT_WG_2D_H
#include "curvedWaveGuide2D.hpp"
#include <complex>
#include <vector>

class ControlFile;

typedef std::complex<double> cdouble;

/** A straight waveguide can be view as a curved one with radius of curvature R = infty */
class StraightWG2D: public CurvedWaveGuideFD
{
public:
  StraightWG2D(): CurvedWaveGuideFD("StraightWaveGuide"){};

  /** Fill JSON object with parameters relevant to this class */
  void fillInfo( Json::Value &obj ) const override final;

  /** Initialize the waveguide. Relevant if loading an old solution */
  void init( const ControlFile &ctl ) override final;

  /** Extract the wavefield at the edge cooresponding to the x-coordinate wcrd */
  void extractField( double wcrd, std::vector<cdouble> &res ) const override final;
protected:
  /** Checks if the coordinates are inside the waveguide */
  bool isInsideGuide( double x, double z ) const override final;

  /** Returns the lower waveguide border at position z */
  double waveGuideStartX( double z ) const override final;

  /** Returns the upper waveguide border at position z */
  double waveGuideEndX( double z ) const override final;
};
#endif
