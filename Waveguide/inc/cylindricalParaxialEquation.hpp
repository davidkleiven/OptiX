#ifndef CYLINDRICAL_PARAXIAL_EQUATION_H
#define CYLINDRICAL_PARAXIAL_EQUATION_H
#include "paraxialEquation.hpp"

class CylindricalParaxialEquation: public ParaxialEquation
{
public:
  CylindricalParaxialEquation(){};
  double F( double r, double theta ) const override final;
  double G( double r, double theta ) const override final;
  double H( double r, double theta ) const override final;
  double J( double r, double theta ) const override final;
  void setRadiusOfCurvature( double Rcurv ){ R = Rcurv; };
  virtual cdouble phaseFactor( double k, double theta ) const override;
private:
    double R{1.0};
};

#endif
