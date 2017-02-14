#ifndef GEOMETRY_SPHERE_H
#define GEOMETRY_SPHERE_H
#include "paraxialSimulation.hpp"

struct Point3D
{
  double x;
  double y;
  double z;
};

class Sphere: public ParaxialSimulation
{
public:
  Sphere( const Point3D &center, double radius );

  /** Return the material properties of the cylinder */
  virtual void getXrayMatProp( double x, double y, double z, double &delta, double &beta ) const override;

  /** Set the cylinder material */
  void setMaterial( const char* name );

  /** Run the simulation without absorption */
  void noAbsorption(){ withAbsorption = false; };
protected:
  Sphere( const Point3D &inCenter, double radius, const char* name );
  Point3D center;
  double radius{1.0};
  double delta{0.0};
  double beta{0.0};
  bool withAbsorption{true};
};

class CoatedSphere: public Sphere
{
public:
  CoatedSphere( const Point3D &center, double radius ): Sphere(center, radius, "CoatedSphere"){};

  /** Set the material of the coating */
  void setCoatingMaterial( const char* name );

  /** Set the thickness of the coating in nano meters */
  void setThickness( double newThickness ){ thickness = newThickness; };

  /** Return the material properties */
  virtual void getXrayMatProp( double x, double y, double z, double &delta, double &beta ) const override;
private:
  double deltaCoat{0.0};
  double betaCoat{0.0};
  double thickness{0.0};
};
#endif
