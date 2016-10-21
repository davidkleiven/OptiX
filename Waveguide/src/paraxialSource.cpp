#include "paraxialSource.hpp"
#include <cmath>

const double PI = acos(-1.0);

void ParaxialSource::setWavelength( double lambda )
{
  wavenumber = 2.0*PI/lambda;
}

double ParaxialSource::getWavelength() const
{
  return 2.0*PI/wavenumber;
}
