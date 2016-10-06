#ifndef NUMEROV_SOLVER_H
#define NUMEROV_SOLVER_H
#include "solver1D.hpp"

struct InitialCondition
{
  double x;
  double value;
};

class Numerov: public Solver1D
{
public:
  Numerov(): Solver1D("Numerov"){};
  void solve() override;
  void setUpperInitialCondition( double x, double value );
  void setLowerInitialCondition( double x, double value );
private:
  double stepsize;
  InitialCondition upper;
  InitialCondition lower;
};
#endif
