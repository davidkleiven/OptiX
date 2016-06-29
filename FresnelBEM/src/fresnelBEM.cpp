#include <cstdlib>
#include <iostream>
#include <libscuff.h>
#include <complex>

int main(int argc, char **argv)
{
  scuff::RWGGeometry geo = scuff::RWGGeometry("square.msh");
  SetLogFileName("fresnel.log");
  geo.SetLogLevel(SCUFF_VERBOSELOGGING);
  
  // Allocate BEM matrix and rhs vector
  HMatrix *matrix = geo.AllocateBEMMatrix();
  HVector *rhsVec = geo.AllocateRHSVector();

  // Source definition
  double sourcePosition[3] = {0.0,0.0,0.3};
  double kHat[3] = {0.0,0.0,-1.0};
  std::complex<double> E0[3] = {1.0,0.0,0.0}; 
  PlaneWave pw(E0, kHat);
  double omega = 1.0;
  
  // Array for storing the fields EH={Ex,Ey,Ez,Hx,Hy,Hz}
  std::complex<double> EH[6];

  delete matrix;
  delete rhsVec;
  return 0;
}
