#include "cladding.hpp"
#include <cmath>

using namespace std;

const double PI = acos(-1.0);
void Cladding::computePotential()
{
  potential = 4.0*PI*thomsonPenetrationDepth*electronDensity;
}

void Cladding::setElectronDensity( double eDensity )
{
  electronDensity = eDensity;
  computePotential();
};
void Cladding::setThomsonPenetrationDepth( double depth )
{
  thomsonPenetrationDepth = depth;
  computePotential();
};
