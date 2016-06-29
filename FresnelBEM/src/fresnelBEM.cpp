#include <cstdlib>
#include <iostream>
#include <libscuff.h>
#include <complex>
#include <cmath>

const double PI = acos(-1.0);
enum class Polarisation_t{S, P};

void getE0_p( const double kHat[3], std::complex<double> E0[3] )
{
  // Assuming the H field is H=(Hx, Hy, Hz) = (1.0, 0.0, 0.0)
  E0[0] = 0.0;
  E0[1] = -kHat[2];
  E0[2] = kHat[1];
}
 
int main(int argc, char **argv)
{
  scuff::RWGGeometry geo = scuff::RWGGeometry("halfspace.scuffgeo");
  SetLogFileName("fresnel.log");
  geo.SetLogLevel(SCUFF_VERBOSELOGGING);
  
  // Allocate BEM matrix and rhs vector
  HMatrix *matrix = geo.AllocateBEMMatrix();
  HVector *rhsVec = geo.AllocateRHSVector();

  // Source definition
  double sourcePosition[3] = {0.0,0.0,0.3};
  double kHat[3] = {0.0,0.0,-1.0};
  std::complex<double> E0_s[3] = {1.0,0.0,0.0}; 
  PlaneWave pw(E0_s, kHat);
  double omega = 1.0;
  
  // Array for storing the fields EH={Ex,Ey,Ez,Hx,Hy,Hz}
  std::complex<double> EH[6];

  double theta = 0.0;
  const double dtheta = 5.0;
  const double thetamax = 4.0;
  Polarisation_t pol[2] = {Polarisation_t::S, Polarisation_t::P};
  double kBloch[3] = {0.0,0.0,0.0};

  std::cout << "Assembling BEM matrix..." << std::endl;
  geo.AssembleBEMMatrix(omega, kBloch, matrix);
  while ( theta < thetamax )
  {
    std::cout << "Theta="<<theta<<std::endl;
    double kz = -cos(theta*PI/180.0);
    double ky = sin(theta*PI/180.0);
    kHat[1] = ky;
    kHat[2] = kz;
    kBloch[1] = ky*omega;
    pw.SetnHat(kHat);

    // Assemble matrix
    matrix->LUFactorize();

    geo.AssembleRHSVector(omega, kBloch, &pw, rhsVec);
    matrix->LUSolve(rhsVec);
    
    theta += dtheta;
  }

  delete matrix;
  delete rhsVec;
  return 0;
}
