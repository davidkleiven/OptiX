#ifndef TRANSPARENT_BOUNDARY_CONDITION_H
#define TRANSPARENT_BOUNDARY_CONDITION_H
#include "boundaryCondition.hpp"

class TransparentBC: public BoundaryCondition
{
public:
  TransparentBC();
  virtual cdouble neighbourCoupling( const Solver2D &solver, double x, double z ) const override;
};
#endif
