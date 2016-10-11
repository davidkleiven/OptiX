#include "curvedWaveGuide2D.hpp"
#include <cassert>
#include <cmath>
#include <vector>
#include "cladding.hpp"
#include "solver2D.hpp"

using namespace std;

cdouble CurvedWaveGuideFD::getRefractiveIndex( double x, double z) const
{
  // Some assertions for debugging
  assert( cladding != NULL );
  assert ( x >= xDisc->min );
  assert ( x <= xDisc->max );
  assert ( z >= zDisc->min );
  assert ( z <= zDisc->max );
  if ( isInsideGuide( x, z) )
  {
    return 1.0;
  }
  return cladding->getRefractiveIndex();
}

bool CurvedWaveGuideFD::isInsideGuide( double x, double z ) const
{
  // TODO: Issues with numerical pressicion here?
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
