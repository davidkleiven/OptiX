#ifndef GAUSSIAN_BEAM_H
#define GAUSSIAN_BEAM_H
#include "paraxialSource.hpp"

class GaussianBeam: public ParaxialSource
{
public:
  GaussianBeam(): ParaxialSource("gaussianBeam"){};
  cdouble get( double x, double z ) const override final;
  void setWaist( double w ){ waist = w; };
  double beamDivergence() const;
  void info( Json::Value &obj ) const override;
private:
  double waist{1.0};

  // Help functions
  double rayleighRange() const;
  double guoyPhase( double z ) const;
  double inverseRadiusOfCurvature( double z ) const;
  double spotSize( double z ) const;
};

#endif
