#ifndef FIXED_VALUES_SOURCE_H
#define FIXED_VALUES_SOURCE_H
#include "paraxialSource.hpp"
#include <armadillo>

/** Class for specifying the sources from known values */
class FixedValuesSource: public ParaxialSource
{
public:
  FixedValuesSource(): ParaxialSource("FixedValuesSource"){};

  /* Evaluates the source at coordinate x, z */
  virtual cdouble get( double x, double z ) const override;

  /** Set limits */
  void setLimits( double xminNew, double xmaxNew ){ xmin=xminNew; xmax=xmaxNew; };

  /** Set the values*/
  void setData( const arma::cx_vec *data ){ values = data; };
private:
  bool limitsSet{false};
  unsigned int indx( double x ) const;
  double getX( unsigned int n ) const;
  double xmin{0.0}, xmax{0.0};
  const arma::cx_vec *values{NULL};
};
#endif
