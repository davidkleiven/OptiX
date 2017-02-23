#ifndef GEOMETRY_SPHERE_H
#define GEOMETRY_SPHERE_H
#include "materialFunction.hpp"
#include <cmath>

struct Point3D
{
  double x;
  double y;
  double z;
};

class Sphere: public MaterialFunction
{
public:
  Sphere( const Point3D &center, double radius );

  /** Return the material properties of the cylinder */
  virtual void getXrayMatProp( double x, double y, double z, double &delta, double &beta ) const override;

  /** Set the cylinder material */
  void setMaterial( const char* name, double energy );

  /** Run the simulation without absorption */
  void noAbsorption(){ withAbsorption = false; };
protected:
  Point3D center;
  double radius{1.0};
  double delta{0.0};
  double beta{0.0};
  bool withAbsorption{true};
};

class CoatedSphere: public Sphere
{
public:
  CoatedSphere( const Point3D &center, double radius, double thickness ): Sphere(center, radius), thickness(thickness){};

  /** Set the material of the coating */
  void setCoatingMaterial( const char* name, double energy );

  /** Return the material properties */
  virtual void getXrayMatProp( double x, double y, double z, double &delta, double &beta ) const override;
private:
  double deltaCoat{0.0};
  double betaCoat{0.0};
  double thickness{0.0};
};
#endif
