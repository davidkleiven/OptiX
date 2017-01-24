#include "cylinder2D.hpp"
#include "refractiveIndex.hpp"
#include <cmath>
#include <iostream>

using namespace std;
Cylinder2D::Cylinder2D( double x0, double z0, double radius ):ParaxialSimulation("Cylinder2D"), \
x0(x0), z0(z0), radius(radius){};

void Cylinder2D::setMaterial( const char* name )
{
  RefractiveIndex refr;
  refr.load(name);
  double energy = getEnergy();
  delta = refr.getDelta( energy );
  beta = refr.getBeta( energy );
}

void Cylinder2D::getXrayMatProp( double x, double z, double &matDelta, double &matBeta ) const
{
  double distanceFromCenterSq = pow( x-x0,2 ) + pow( z-z0, 2 );
  if ( distanceFromCenterSq < radius*radius )
  {
    matDelta = delta;
    matBeta = beta;
    return;
  }
  matDelta = 0.0;
  matBeta = 0.0;
}
