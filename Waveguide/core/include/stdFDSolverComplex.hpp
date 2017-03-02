#ifndef STANDARD_FD_SOLVER_COMPLEX_H
#define STANDARD_FD_SOLVER_COMPLEX_H

/** Class for handling standard central differences, but wita complex values solution */
class StdFDSolverComplex: public Solver1D
{
public:
  StdFDSolverComplex(): Solver1D("StandardFDComplex"){};

  /** Solve the equation */
  void solve() override final;

  /** Set the stepsize */
  void setStepsize( double newstep ){ stepsize = newstep; };

  /** Fill JSON object with parameters relevant to this class */
  void fillJsonObj( Json::Value &obj ) const override final;
private:
  double stepsize;
};
