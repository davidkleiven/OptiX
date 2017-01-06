#include "transparentBoundaryCondition.hpp"
#include "solver2D.hpp"

TransparentBC::TransparentBC():BoundaryCondition("TransparentBC"){};

cdouble TransparentBC::neighbourCoupling( const Solver2D &solver, double x, double z ) const
{
  double xmin = solver.getSimulator().transverseDiscretization().min;
  double xmax = solver.getSimulator().transverseDiscretization().max;

  cdouble ratio = 1.0;
  // Check which boundary is closer
  if ( abs(x-xmin) < abs(x-xmax) )
  {
    // Lower boundary --> smallest induces
    ratio = solver->getLastSolution()(0)/solver->getLastSolution(1);
  }
  else
  {
    // Upper boundary --> largest indices
    unsigned int length = solver->getLastSolution().n_elem;
    ratio = solver->getLastSolution()(length-1)/solver->getLastSolution(length-2);
  }

  cdouble im(0.0,1.0);
  cdouble kdx = log(ratio)/im;

  // Check that the real part of k is positive
  if ( kdx.real() < 0.0 ) kdx.real(0.0);

  return exp(im*kdx);
}
