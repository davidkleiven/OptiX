#include "numerovSolver.hpp"

void Numerov::setUpperInitialCondition( double x, double value )
{
  upper.x = x;
  upper.value = value;
}

void Numerov::setLowerInitialCondition( double x, double value )
{
  lower.x = x;
  lower.value = value;
}
