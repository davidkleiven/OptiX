#ifndef CYLINDER_2D_H
#define CYLINDER_2D_H
#include <complex>
#include "paraxialSimulation.hpp"

typedef std::complex<double> cdouble;

/** Class for solving the paraxial waveequation for a cylinder */
class Cylinder2D: public ParaxialSimulation
{
public:
  Cylinder2D( double x0, double z0, double radius );

  /** Return the material properties of the cylinder */
  void getXrayMatProp( double x, double z, double &delta, double &beta ) const override;

  /** Set the cylinder material */
  void setMaterial( const char* name );

  /** Set the transverse boundary conditions. Here this is simple a planewave */
  cdouble transverseBC( double z ) const override { return 1.0; };

  /** Override the transverse boundary conditions */
  cdouble transverseBC( double z, Boundary_t bnd ) const override { return 1.0; };
private:
  double x0{0.0};
  double z0{0.0};
  double radius{0.0};
  double delta{0.0};
  double beta{0.0};
};
#endif