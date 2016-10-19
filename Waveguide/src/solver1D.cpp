#include "solver1D.hpp"
#include <cassert>

Solver1D::~Solver1D()
{
  delete solution;
}

void Solver1D::fillJsonObj( Json::Value &obj ) const
{
  obj["name"] = name;
}

unsigned int Solver1D::closestIndx( double x ) const
{
  assert ( x >= x1 );
  assert ( x < x2 );
  unsigned int N = getEigenVectorSize();
  return N*( x- x1 )/( x2-x1 );
}

double Solver1D::getSolution( double x, unsigned int eigenmode ) const
{
  if (( x < x1) || ( x > x2 ))
  {
    return 0.0;
  }
  unsigned int indx = closestIndx(x);
  double dx = (x2-x1)/static_cast<double>(getEigenVectorSize());
  double x_low = x1+dx*indx;
  double x_high = x1+dx;
  return ((x-x_low)*getSolution()(indx,eigenmode) + (x_high-x)*getSolution()(indx+1,eigenmode))/dx;
}
