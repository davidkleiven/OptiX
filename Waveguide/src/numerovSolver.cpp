#include "numerovSolver.hpp"
#include <gsl/gsl_errno.h>
#include <cassert>
#include <iostream>
#include "waveGuide.hpp"
#include "cladding.hpp"
#include <stdexcept>
using namespace std;

void Numerov::setUpperInitialCondition( double x1, double value1, double value2 )
{
  upper.x1 = x1;
  upper.value1 = value1;
  upper.value2 = value2;

}

void Numerov::setLowerInitialCondition( double x1, double value1, double value2 )
{
  lower.x1 = x1;
  lower.value1 = value1;
  lower.value2 = value2;

}

void Numerov::setPropgationWavenumberLimits( double beta1, double beta2 )
{
  if ( waveguide == NULL )
  {
    throw (runtime_error("No waveguide set\n"));
  }
  beta_min = beta1;//*waveguide->getWavenumber();
  beta_max = beta2;//*waveguide->getWavenumber();
  initPropagationWavenumberLimits = false;
}

void Numerov::iterateForward( unsigned int last )
{
  assert( last >= 1 );
  double currX = lower.x1 + static_cast<double>(last)*stepsize;
  double prevX = currX - stepsize;
  double nextX = currX + stepsize;

  double yn = (*solution)[last];
  double yn_m1 = (*solution)[last-1];
  (*solution)[last + 1] = 2.0*yn*alpha_n(currX) - yn_m1*alpha_nm1(prevX);
  (*solution)[last + 1] /= alpha_np1(nextX);

  /*
  if ( (*solution)[last] < 0.0 )
  {
    cout << "x:" << currX << " value:" << (*solution)[last] << endl;
  }
  */
}

void Numerov::iterateBackward( unsigned int last )
{
  assert( last <= solution->size()-2);
  double currX = lower.x1 + static_cast<double>(last)*stepsize;
  double prevX = currX + stepsize;
  double nextX = currX - stepsize;

  double yn = (*solution)[last];
  double yn_p1 = (*solution)[last+1];
  (*solution)[last-1] = 2.0*yn*alpha_n(currX) - yn_p1*alpha_np1(prevX);
  (*solution)[last-1] /= alpha_nm1(nextX);
}

double Numerov::effectivePotential( double x ) const
{
  double k = waveguide->getWavenumber();
  //return -(waveguide->potential(x) + eigenvalue*eigenvalue - k*k);
  //cout << "Potential(:" << waveguide->potential(x) - eigenvalue << ") ";
  return -(waveguide->potential(x) - currentEigval);
}

double Numerov::alpha_np1( double x ) const
{
  return 1.0 + stepsize*stepsize*effectivePotential(x)/12.0;
}

double Numerov::alpha_n( double x ) const
{
  return (1.0 - 5.0*stepsize*stepsize*effectivePotential(x)/12.0);
}

double Numerov::alpha_nm1( double x ) const
{
  return 1.0 + stepsize*stepsize*effectivePotential(x)/12.0;
}

void Numerov::iterateAll()
{
  unsigned int N = solution->size();
  unsigned int meetingPoint = N/2;
  for ( unsigned int i=1;i<meetingPoint-1;i++ )
  {
    iterateForward(i);
  }
  for ( unsigned int i=solution->size()-2;i>meetingPoint;i--)
  {
    iterateBackward(i);
  }
}

double Numerov::rootSolverFunction( double beta, void *params )
{
  Numerov* self = static_cast<Numerov*>(params);
  self->currentEigval = beta;
  self->iterateAll();
  unsigned int N = self->solution->size();
  unsigned int middle = N/2;
  return (*self->solution)[middle-1] - (*self->solution)[middle];
}

void Numerov::solve()
{
  if ( stepsize < 0.0 )
  {
    throw (runtime_error("Invalid stepsize"));
  }

  cerr << beta_min << " " << beta_max << endl;
  unsigned int N = (upper.x1 + stepsize - lower.x1)/stepsize;

  // Set the solution size
  solution->clear();
  solution->resize(N);
  fill( solution->begin(), solution->end(), 0.0);
  (*solution)[0] = lower.value1;
  (*solution)[1] = lower.value2;
  (*solution)[solution->size()-2] = upper.value1;
  (*solution)[solution->size()-1] = upper.value2;

  gsl_function F;
  F.function = &rootSolverFunction;
  F.params = this;

  if ( initPropagationWavenumberLimits )
  {
    beta_min = 0.0;
    beta_max = 0.99*waveguide->getCladding().getPotential();
  }

  const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

  gsl_root_fsolver_set(s, &F, beta_min, beta_max);

  int status = GSL_CONTINUE;
  iter = 0;
  while (( status == GSL_CONTINUE ) && ( iter < maxIterations ))
  {
    //cout << "\n=================== ITER " << iter << " =====================================\n";
    iter++;
    currentEigval = gsl_root_fsolver_root(s);
    beta_min = gsl_root_fsolver_x_lower(s);
    beta_max = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(beta_min, beta_max, 0.0, 0.00001);
  }
  gsl_root_fsolver_free(s);
  if ( iter == maxIterations )
  {
    cout << "Warning: Reached max number of iterations...\n";
  }
}

void Numerov::fillJsonObj( Json::Value &obj ) const
{
  Solver1D::fillJsonObj(obj);
  obj["xmin"] = lower.x1;
  obj["xmax"] = upper.x1+stepsize;
  obj["iterations"] = iter;
  obj["maxIterations"] = maxIterations;
  obj["beta_min"] = beta_min;
  obj["beta_max"] = beta_max;
}
