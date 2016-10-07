#include "numerovSolver.hpp"
#include <gsl/gsl_errno.h>
#include <cassert>
#include <iostream>
#include "waveGuide.hpp"
using namespace std;

void Numerov::setUpperInitialCondition( double x1, double value1, double x2, double value2 )
{
  upper.x1 = x1;
  upper.value1 = value1;
  upper.x2 = x2;
  upper.value2 = value2;

}

void Numerov::setLowerInitialCondition( double x1, double value1, double x2, double value2 )
{
  lower.x1 = x1;
  lower.value1 = value1;
  lower.x2 = x2;
  lower.value2 = value2;

}

void Numerov::setPropgationWavenumberLimits( double beta1, double beta2 )
{
  beta_min = beta1;
  beta_max = beta2;
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
  (*solution)[last+1] /= alpha_np1(nextX);
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
  return waveguide->potential(x) - eigenvalue*eigenvalue + k*k;
}

double Numerov::alpha_np1( double x ) const
{
  return 1.0 + stepsize*stepsize*effectivePotential(x)/12.0;
}

double Numerov::alpha_n( double x ) const
{
  return 2.0*(1.0 - 5.0*stepsize*stepsize*effectivePotential(x)/12.0);
}

double Numerov::alpha_nm1( double x ) const
{
  return 1.0 + stepsize*stepsize*effectivePotential(x)/12.0;
}

void Numerov::iterateAll()
{
  unsigned int N = solution->size();
  for ( unsigned int i=1;i<N/2;i++ )
  {
    iterateForward(i);
  }
  for ( unsigned int i=solution->size()-2;i>=N/2;i--)
  {
    iterateBackward(i);
  }
}

double Numerov::rootSolverFunction( double beta, void *params )
{
  Numerov* self = static_cast<Numerov*>(params);
  self->eigenvalue = beta;
  self->iterateAll();
  unsigned int N = self->solution->size();
  unsigned int middle = N/2;
  return (*self->solution)[middle-1] - (*self->solution)[middle];
}

void Numerov::solve()
{
  double middle = 0.5*(upper.x2 + lower.x1);
  unsigned int N = (upper.x2 - lower.x1)/stepsize;

  // Set the solution size
  solution->clear();
  solution->resize(N);
  gsl_function F;
  F.function = &rootSolverFunction;
  F.params = this;

  if ( initPropagationWavenumberLimits )
  {
    beta_min = 0.8*waveguide->getWavenumber();
    beta_max = 1.0*waveguide->getWavenumber();
  }

  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

  gsl_root_fsolver_set(s, &F, beta_min, beta_max);

  int status = GSL_CONTINUE;
  iter = 0;
  while (( status == GSL_CONTINUE ) && ( iter < maxIterations ))
  {
    iter++;
    eigenvalue = gsl_root_fsolver_root(s);
    beta_min = gsl_root_fsolver_x_lower(s);
    beta_max = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(beta_min, beta_max, 0.0, 0.001);
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
  obj["xmax"] = upper.x2;
  obj["iterations"] = iter;
  obj["maxIterations"] = maxIterations;
  obj["beta_min"] = beta_min;
  obj["beta_max"] = beta_max;
}
