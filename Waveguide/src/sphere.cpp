#include "sphere.hpp"
#include "refractiveIndex.hpp"
#include <iostream>

using namespace std;

Sphere::Sphere( const Point3D &inCenter, double radius ): Sphere(inCenter, radius, "SingleSpehre" ){};

Sphere::Sphere( const Point3D &inCenter, double radius, const char* name ): ParaxialSimulation(name), radius(radius)
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

void Sphere::setGroupAttributes()
{
  if ( maingroup == NULL ) return;

  ParaxialSimulation::setGroupAttributes();

  H5::DataSpace attribSpace(H5S_SCALAR);

  H5::Attribute att = maingroup->createAttribute( "delta", H5::PredType::NATIVE_DOUBLE, attribSpace );
  att.write( H5::PredType::NATIVE_DOUBLE, &delta );

  att = maingroup->createAttribute( "beta", H5::PredType::NATIVE_DOUBLE, attribSpace );
  att.write( H5::PredType::NATIVE_DOUBLE, &beta );

  att = maingroup->createAttribute( "radius", H5::PredType::NATIVE_DOUBLE, attribSpace );
  att.write( H5::PredType::NATIVE_DOUBLE, &radius );

  att = maingroup->createAttribute( "withAbsorption", H5::PredType::NATIVE_HBOOL, attribSpace );
  att.write( H5::PredType::NATIVE_HBOOL, &withAbsorption );
}

////////////////////////////////////////////////////////////////////////////////
void CoatedSphere::getXrayMatProp( double x, double y, double z, double &matDelta, double &matBeta ) const
{
  Sphere::getXrayMatProp( x, y, z, matDelta, matBeta );
  const double ZERO = 1E-12;
  if (( matDelta > ZERO ) || ( matBeta > ZERO ))
  {
    // Inside the core of the sphere
    return;
  }

  double distSq = pow( x - center.x, 2 ) + pow( y - center.y, 2 ) + pow( z - center.z, 2 );

  if (( distSq < pow( radius+thickness, 2) ) && ( distSq > pow(radius,2) ))
  {
    matDelta = deltaCoat;
    matBeta = betaCoat;
    return;
  }

  matDelta = 0.0;
  matBeta = 0.0;
}

void CoatedSphere::setCoatingMaterial( const char* name )
{
  RefractiveIndex refr;
  refr.load(name);
  double energy = getEnergy();
  deltaCoat = refr.getDelta( energy );

  if ( withAbsorption )
  {
    betaCoat = refr.getBeta( energy );
  }
  else
  {
    betaCoat = 0.0;
  }
}

void CoatedSphere::setGroupAttributes()
{
  if ( maingroup == NULL ) return;
  Sphere::setGroupAttributes();

  H5::DataSpace attribSpace(H5S_SCALAR);

  H5::Attribute att = maingroup->createAttribute( "deltaCoating", H5::PredType::NATIVE_DOUBLE, attribSpace );
  att.write( H5::PredType::NATIVE_DOUBLE, &deltaCoat );

  att = maingroup->createAttribute( "betaCoating", H5::PredType::NATIVE_DOUBLE, attribSpace );
  att.write( H5::PredType::NATIVE_DOUBLE, &betaCoat );

  att = maingroup->createAttribute( "coatThickness", H5::PredType::NATIVE_DOUBLE, attribSpace );
  att.write( H5::PredType::NATIVE_DOUBLE, &thickness );
}
