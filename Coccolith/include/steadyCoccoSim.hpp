#ifndef STEADY_COCCOLITH_SIM_H
#define STEADY_COCCOLITH_SIM_H
#include "coccolithSim.hpp"

class SteadyCoccolithSim: public CoccolithSimulation
{
public:
  SteadyCoccolithSim(){};

  /** Initialize the source */
  void initSource( double freq );

  /** Runs the simulation */
  virtual void run() override;

  /** Run in time domain for a given time */
  void stepForTime( double sec );

  /** Export results */
  virtual void exportResults() override;

  /** Initializes the simulation */
  virtual void init() override;


  /** Tolerance for convergence of the iterative solver */
  double tolerance{1E-8}; // MEEP default value

  /** L-parameter in BiCStab-L algorithm */
  int BiCStabL{2};               // MEEP default value

  /** Maximum number of iterations */
  int maxiters{10000};    // MEEP default value
private:
  meep::continuous_src_time *contSource{NULL};
};
#endif
