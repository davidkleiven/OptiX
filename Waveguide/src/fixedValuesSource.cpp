#include "fixedValuesSource.hpp"
#include <cassert>

cdouble FixedValuesSource::get( double x, double z ) const
{
  assert( values != NULL );
  if ( x < xmin ) return (*values)(0);
  else if ( x >= xmax ) return (*values)(values->n_elem-1);

  unsigned int n = indx( x );
  if ( n >= values->n_elem-1 ) return (*values)(values->n_elem-1);

  double xlow = getX(n);
  double xhigh = getX(n+1);
  double w = (x-xlow)/(xhigh-xlow);
  return (*values)(n+1)*w + (1.0-w)*(*values)(n);
}

unsigned int FixedValuesSource::indx( double x ) const
{
  return (x-xmin)*values->n_elem/(xmax-xmin);
}

double FixedValuesSource::getX( unsigned int n ) const
{
  double dx = (xmax-xmin)/values->n_elem;
  return n*dx+xmin;
}
