#include "waveGuideFDSimulation.hpp"
#include "solver2D.hpp"
#include "cladding.hpp"

WaveGuideFDSimulation::WaveGuideFDSimulation(): xDisc(new Disctretization), zDisc(new Disctretization), name(""){};

WaveGuideFDSimulation::WaveGuideFDSimulation( const char* wgname ): WaveGuideFDSimulation()
{
  name = wgname;
}

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
  return (xDisc->max - xDisc->min)/xDisc->step + 1.0;
}

unsigned int WaveGuideFDSimulation::nodeNumberLongitudinal() const
{
  return ( zDisc->max - zDisc->min)/zDisc->step + 1.0;
}

void WaveGuideFDSimulation::setSolver( Solver2D &solv )
{
  solver = &solv;
  solver->setGuide( *this );
}

void WaveGuideFDSimulation::setCladding( const Cladding &clad )
{
  cladding = &clad;
}
