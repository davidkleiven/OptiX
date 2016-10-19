#include "linearMap1D.hpp"
#include <iostream>

using namespace std;

void LinearMap1D::initialize( const CorrespondingPoints &p1, const CorrespondingPoints &p2 )
{
  jacobian = (p2.xFrom - p1.xFrom)/(p2.xTo - p1.xTo);
  shift = p1.xTo - jacobian*p1.xFrom;
}

double LinearMap1D::get( double xFrom ) const
{
  return jacobian*xFrom + shift;
}
