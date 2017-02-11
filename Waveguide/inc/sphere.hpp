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
  void getXrayMatProp( double x, double y, double z, double &delta, double &beta ) const override;

  /** Set the cylinder material */
  void setMaterial( const char* name );

  /** Run the simulation without absorption */
  void noAbsorption(){ withAbsorption = false; };
private:
  Point3D center;
  double radius{1.0};
  double delta{0.0};
  double beta{0.0};
  bool withAbsorption{true};
};
#endif
