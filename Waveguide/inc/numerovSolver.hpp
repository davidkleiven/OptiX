#ifndef NUMEROV_SOLVER_H
#define NUMEROV_SOLVER_H
#include "solver1D.hpp"
#include <gsl/gsl_roots.h>

struct InitialCondition
{
  double x1;
  double value1;
  double value2; // Value at x1+stepsize
};

class Numerov: public Solver1D
{
public:
  Numerov(): Solver1D("Numerov"){};
  void solve() override final;
  void setUpperInitialCondition( double x1, double value1, double value2 );
  void setLowerInitialCondition( double x1, double value1, double value2 );
  void setInitialEigenvalue( double eigen ){ eigenvalue = eigen; };
  void setMaxIter( unsigned int mxIter ){ maxIterations = mxIter; };
  void setPropgationWavenumberLimits( double beta_min, double beta_max );
  void fillJsonObj( Json::Value &obj ) const override final;
  void setStepsize( double step ){ stepsize=step; };
protected:
  double stepsize{-1.0};
  unsigned int iter{0};
  InitialCondition upper;
  InitialCondition lower;
  void iterateForward( unsigned int last );
  void iterateBackward( unsigned int last );
  double beta_min, beta_max; // Limits for the propagation wave number
  void iterateAll();
  unsigned int maxIterations{100};
  bool initPropagationWavenumberLimits{true};

  // Helpers for computing different terms of the Numerov method
  double alpha_np1( double x ) const;
  double alpha_n( double x ) const;
  double alpha_nm1( double x ) const;
  double effectivePotential( double x ) const;

  static double rootSolverFunction( double beta, void *params ); // GSL need this
};
#endif
