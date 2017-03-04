#ifndef STEADY_COCCOLITH_SIM_H
#define STEADY_COCCOLITH_SIM_H
#include "coccolithSim.hpp"

class SteadyCoccolithSim: public CoccolithSimulation
{
public:
  SteadyCoccolithSim(){};
  virtual ~SteadyCoccolithSim();

  /** Initialize the source */
  void initSource( double freq );

  /** Runs the simulation */
  virtual void run() override;

  /** Export results */
  virtual void exportResults() override;
private:
  meep::continuous_src_time *contSource{NULL};
};
#endif
