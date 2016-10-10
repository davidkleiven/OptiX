#include "waveGuideFDSimulation.hpp"

WaveGuideFDSimulation::WaveGuideFDSimulation(): xDisc(new Disctretization), zDisc(new Disctretization){};

WaveGuideFDSimulation::~WaveGuideFDSimulation()
{
  delete xDisc;
  delete zDisc;
}
void WaveGuideFDSimulation::setTransverseDiscretization( double xmin, double xmax, double step )
{
  xDisc->min = xmin;
  xDisc->max = xmax;
  xDisc->step = step;
}

void WaveGuideFDSimulation::setLongitudinalDiscretization( double xmin, double xmax, double step )
{
  zDisc->min = xmin;
  zDisc->max = xmax;
  zDisc->step = step;
}

unsigned int WaveGuideFDSimulation::nodeNumberTransverse() const
{
  return (xDisc->max - xDisc->min)/xDisc->step;
}

unsigned int WaveGuideFDSimulation::nodeNumberLongitudinal() const
{
  return ( zDisc->max - zDisc->min)/zDisc->step;
}
