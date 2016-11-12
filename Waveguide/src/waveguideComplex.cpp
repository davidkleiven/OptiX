#include "waveguideComplex.hpp"

WaveguideComplex::WaveguideComplex():WaveGuideLargeCurvature("LargeRadiusOfCurvatureComplex")
{
  useComplexPotential = true;
}

cdouble WaveguideComplex::complexPotential( double x ) const
{
  double beta = cladding->getBeta();
  double delta = cladding->getDelta();
}
