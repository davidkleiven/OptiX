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
};

#endif
