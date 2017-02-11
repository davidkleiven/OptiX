#include "sphere.hpp"
#include "refractiveIndex.hpp"
#include <iostream>

using namespace std;

Sphere::Sphere( const Point3D &inCenter, double radius ): ParaxialSimulation("SingleSphere"), radius(radius)
{
  center.x = inCenter.x;
  center.y = inCenter.y;
  center.z = inCenter.z;
}

void Sphere::setMaterial( const char* name )
{
  RefractiveIndex refr;
  refr.load(name);
  double energy = getEnergy();
  delta = refr.getDelta( energy );

  if ( withAbsorption )
  {
    beta = refr.getBeta( energy );
  }
  else
  {
    beta = 0.0;
  }
}

void Sphere::getXrayMatProp( double x, double y, double z, double &matDelta, double &matBeta ) const
{
  double distSq = pow( x - center.x, 2 ) + pow( y - center.y, 2 ) + pow( z - center.z, 2 );

  if ( distSq < radius*radius )
  {
    matDelta = delta;
    matBeta = beta;
    return;
  }

  matDelta = 0.0;
  matBeta = 0.0;
}
