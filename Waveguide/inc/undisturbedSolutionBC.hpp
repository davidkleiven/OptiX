#ifndef UNDISTURBED_SOLUTION_BC_H
#define UNDISTURBED_SOLUTION_BC_H
#include "boundaryCondition.hpp"
#include<complex>

typedef std::complex<double> cdouble;
class ParaxialSource;
class Solver2D;

class UndisturbedSolutionBC: public BoundaryCondition
{
public:
  UndisturbedSolutionBC( const ParaxialSource &src, double delta, double beta );
  virtual cdouble fixedField( double x, double z ) const override;
private:
  const ParaxialSource *src;
  double delta, beta;
};
#endif
