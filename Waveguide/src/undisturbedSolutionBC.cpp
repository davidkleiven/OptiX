#include "paraxialSource.hpp"
#include "undisturbedSolutionBC.hpp"
#include <cmath>

UndisturbedSolutionBC::UndisturbedSolutionBC( const ParaxialSource &src, double delta, double beta ):
BoundaryCondition("UndisturbedBC"), src(&src), delta(delta), beta(beta){};

cdouble UndisturbedSolutionBC::fixedField( double x, double z ) const
{
  cdouble im(0.0,1.0);
  cdouble amplitude = src->get(x,0.0);
  double k = src->getWavenumber();
  return exp(-im*k*delta)*exp(-k*beta*z);
}
