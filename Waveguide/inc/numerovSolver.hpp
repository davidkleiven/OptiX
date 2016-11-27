#ifndef NUMEROV_SOLVER_H
#define NUMEROV_SOLVER_H
#include "solver1D.hpp"
#include <gsl/gsl_roots.h>

/** Struct for setting the initial conditions */
struct InitialCondition
{
  double x1;
  double value1;
  double value2; // Value at x1+stepsize
};

/** Class for implementing a Numerov solver to the 1D Schrodinger equation. Depricated */
class Numerov: public Solver1D
{
public:
  Numerov(): Solver1D("Numerov"){};

  /** Solve the system with the shooting method */
  void solve() override final;

  /** Set initial conditions at x=xmax*/
  void setUpperInitialCondition( double x1, double value1, double value2 );

  /** Set initial conditions at x=xmin */
  void setLowerInitialCondition( double x1, double value1, double value2 );

  /** Set initial guess for the eigenvalue */
  void setInitialEigenvalue( double eigen ){ currentEigval = eigen; };

  /** Set maximum number of iterations */
  void setMaxIter( unsigned int mxIter ){ maxIterations = mxIter; };

  /** Set the domain over valid propagation constants */
  void setPropgationWavenumberLimits( double beta_min, double beta_max );

  /** Fill JSON object with parameters relevant to this class */
  void fillJsonObj( Json::Value &obj ) const override final;

  /** Set the stepsize in nano meters */
  void setStepsize( double step ){ stepsize=step; };
protected:
  double currentEigval;
  double stepsize{-1.0};
  unsigned int iter{0};
  InitialCondition upper;
  InitialCondition lower;

  /** Perform forward iteration */
  void iterateForward( unsigned int last );

  /** Perform backward iteration */
  void iterateBackward( unsigned int last );
  double beta_min, beta_max; // Limits for the propagation wave number

  /** Perform iteration both ways */
  void iterateAll();
  unsigned int maxIterations{100};
  bool initPropagationWavenumberLimits{true};

  // Helpers for computing different terms of the Numerov method
  /** Return coefficient in front of psi_{n+1} */
  double alpha_np1( double x ) const;

  /** Return coefficient in front of psi_n */
  double alpha_n( double x ) const;

  /** Return coefficient in front of psi_{n-1}*/
  double alpha_nm1( double x ) const;

  /** Return the effective potantial */
  double effectivePotential( double x ) const;

  /** Function that should have root when the correct solution is found. Required by GSL*/
  static double rootSolverFunction( double beta, void *params ); // GSL need this
};
#endif
