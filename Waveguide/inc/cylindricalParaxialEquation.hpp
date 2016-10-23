#ifndef CYLINDRICAL_PARAXIAL_EQUATION_H
#define CYLINDRICAL_PARAXIAL_EQUATION_H
#include "paraxialEquation.hpp"

class CylindricalParaxialEquation: public ParaxialEquation
{
public:
  CylindricalParaxialEquation(){};
  cdouble F( double r, double theta ) const override final;
  cdouble G( double r, double theta ) const override final;
  cdouble H( double r, double theta ) const override final;
};

#endif
