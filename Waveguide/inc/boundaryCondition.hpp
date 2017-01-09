#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H
#include <complex>
#include <string>

class Solver2D;

typedef std::complex<double> cdouble;

/** Interface class for transverse boundary conditions */
class BoundaryCondition
{
public:
  BoundaryCondition( const char* name ):name(name){};

  /** Returns the field at the boundary */
  virtual cdouble fixedField( double x, double z ) const{ return 0.0; };

  /** Returns an additional weight that will be added to the neighbouring node in the FD stencil */
  virtual cdouble neighbourCoupling( const Solver2D &solver, double x, double z ) const { return 0.0; };
private:
  std::string name;
};
#endif
