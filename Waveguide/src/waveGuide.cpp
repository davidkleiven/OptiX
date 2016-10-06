#include "waveGuide.hpp"
#include "cladding.hpp"

void WaveGuide1DSimulation::setCladding( const Cladding &newCladding )
{
  cladding = &newCladding;
}
