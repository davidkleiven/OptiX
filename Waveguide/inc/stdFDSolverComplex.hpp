#ifndef STANDARD_FD_SOLVER_COMPLEX_H
#define STANDARD_FD_SOLVER_COMPLEX_H

class StdFDSolverComplex: public Solver1D
{
public:
  StdFDSolverComplex(): Solver1D("StandardFDComplex"){};
  void solve() override final;
  void setStepsize( double newstep ){ stepsize = newstep; };
  void fillJsonObj( Json::Value &obj ) const override final;
private:
  double stepsize;
};
