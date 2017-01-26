#ifndef GAUSSIAN_BEAM_H
#define GAUSSIAN_BEAM_H
#include "paraxialSource.hpp"

/** Gaussian beam source */
class GaussianBeam: public ParaxialSource
{
public:
  GaussianBeam(): ParaxialSource("gaussianBeam"){};

  /** Evaluate the beam at position x, z. Both in nano meters */
  cdouble get( double x, double z ) const override final;

  /** Gaussian beam in 3D */
  cdouble get( double x, double y, double z ) const override final;

  /** Sets the waist in nano meter */
  void setWaist( double w ){ waist = w; };

  /** Get the beam divergence */
  double beamDivergence() const;

  /** Fills a JSON object with parameters specific to this class */
  void info( Json::Value &obj ) const override;
private:
  double waist{1.0};

  // Help functions
  /** Get the rayleigh range */
  double rayleighRange() const;

  /** Get the Guoy phase shift */
  double guoyPhase( double z ) const;

  /** Get the inverse radius of curvature */
  double inverseRadiusOfCurvature( double z ) const;

  /** Get the spot size */
  double spotSize( double z ) const;
};

#endif
