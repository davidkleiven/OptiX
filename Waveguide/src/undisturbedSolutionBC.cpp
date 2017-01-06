#include "paraxialSource.hpp"
#include "undisturbedSolutionBC.hpp"
#include <cmath>

UndisturbedSolutionBC::UndisturbedSolutionBC( const ParaxialSource &src, double delta, double beta ):
BoundaryCondition("UndisturbedBC"), src(&src), delta(delta), beta(beta){};

cdouble UndisturbedSolutionBC::fixedField( const Solver2D &solver, double x, double z ) const
{
  cdouble amplitude = src->get(x,0.0);
  double k = src->getWavenumber();
  return exp(-im*k*delta)*exp(-k*beta*z);
}
