#include "gaussianWG.hpp"
#include "cladding.hpp"
#include <cmath>

using namespace std;

void GaussianWG::getXrayMatProp( double x, double z, double &delta, double &beta ) const
{
  delta = cladding->getDelta()*profile(x,z);
  beta = cladding->getBeta()*profile(x,z);
}

double GaussianWG::profile( double x, double z ) const
{
  double rSq = pow( x+0.5*z*z/R - 0.5*width, 2);
  return 1.0 - exp(-rSq/(2.0*width*width));
}
