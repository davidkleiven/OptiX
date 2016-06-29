#include <cstdlib>
#include <iostream>
#include <libscuff.h>

int main(int argc, char **argv)
{
  scuff::RWGGeometry geo = scuff::RWGGeometry("square.msh");
  SetLogFileName("fresnel.log");
  geo.SetLogLevel(SCUFF_VERBOSELOGGING);
  
  // Allocate BEM matrix and rhs vector
  HMatrix *matrix = geo.AllocateBEMMatrix();
  HVector *rhsVec = geo.AllocateRHSVector();

  double sourcePosition[3] = {0.0,0.0,0.3};

  delete matrix;
  delete rhsVec;
  return 0;
}
