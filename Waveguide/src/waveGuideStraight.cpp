#include "waveGuideStraight.hpp"
#include "cladding.hpp"
#include <iostream>

using namespace std;

double WaveGuideStraight::potential( double x ) const
{
  if (( x > 0.0 ) || ( x < -width ))
  {
    return cladding->getPotential();
  }
  return 0.0;
}

WaveGuideStraight::~WaveGuideStraight(){};
