#ifndef STANDARD_FINITE_DIFFERENCE_SOLVER_H
#define STANDARD_FINITE_DIFFERENCE_SOLVER_H
#include "solver1D.hpp"

class StandardFD: public Solver1D
{
public:
  StandardFD(): Solver1D("StandardFD"){};
  void solve() override final;
  void setLimits( double x0, double x2 );
  void setStepsize( double newstep ){ stepsize = newstep; };
  void fillJsonObj( Json::Value &obj ) const override final;
private:
  double stepsize;
};
#endif
