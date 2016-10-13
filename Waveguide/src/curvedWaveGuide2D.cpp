#include "curvedWaveGuide2D.hpp"
#include <cassert>
#include <cmath>
#include <vector>
#include "cladding.hpp"
#include "solver2D.hpp"
#include <cmath>

using namespace std;

void CurvedWaveGuideFD::getXrayMatProp( double x, double z, double &delta, double &beta) const
{
  // Some assertions for debugging
  assert( cladding != NULL );
  assert ( x >= xDisc->min );
  assert ( x <= xDisc->max );
  assert ( z >= zDisc->min );
  assert ( z <= zDisc->max );
  if ( isInsideGuide( x, z) )
  {
    beta = 0.0;
    delta = 0.0;
    return;
  }
  delta = cladding->getDelta();
  beta = cladding->getBeta();
}

bool CurvedWaveGuideFD::isInsideGuide( double x, double z ) const
{
//  double d = sqrt(x*x + z*z);
  return (2.0*x*R+z*z > 0.0 ) && (2.0*x*R+z*z < 2.0*width*R);
  //return ( d < R+width) && (d > R );
}

void CurvedWaveGuideFD::setBoundaryConditions()
{
  unsigned int Nx = nodeNumberTransverse();
  vector<cdouble> values(Nx, 1.0);
  solver->setLeftBC(&values[0]);
}

void CurvedWaveGuideFD::fillInfo( Json::Value &obj ) const
{
  obj["RadiusOfCurvature"] = R;
  obj["Width"] = width;
}

cdouble CurvedWaveGuideFD::transverseBC( double z ) const
{
  double delta = cladding->getDelta();
  double beta = cladding->getBeta();
  cdouble im(0.0,1.0);
  return exp(-beta*wavenumber*z)*exp(-im*delta*wavenumber*z);
}

void CurvedWaveGuideFD::computeTransmission( double step ) const
{
  // Integrate across waveguide
  double fluxAtZero = 0.0;
  unsigned int wgStart = 0.0;
}
