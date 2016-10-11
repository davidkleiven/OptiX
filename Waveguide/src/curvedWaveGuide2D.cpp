#include "curvedWaveGuide2D.hpp"
#include <cassert>
#include <cmath>
#include <vector>
#include "cladding.hpp"
#include "solver2D.hpp"

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
  }
  delta = cladding->getDelta();
  beta = cladding->getBeta();
}

bool CurvedWaveGuideFD::isInsideGuide( double x, double z ) const
{
  double d = sqrt(x*x + z*z);
  return ( d < R+width) && (d > R );
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
