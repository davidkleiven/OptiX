#ifndef STANDARD_FINITE_DIFFERENCE_SOLVER_H
#define STANDARD_FINITE_DIFFERENCE_SOLVER_H
#include "solver1D.hpp"

/** Class for handling standard central differences solution */
class StandardFD: public Solver1D
{
public:
  StandardFD(): Solver1D("StandardFD"){};

  /** Solve the equation */
  void solve() override final;

  /** Set the stepsize */
  void setStepsize( double newstep ){ stepsize = newstep; };

  /** Fill JSON object with parameters specific to this equation */
  void fillJsonObj( Json::Value &obj ) const override final;
private:
  double stepsize;
};
#endif
