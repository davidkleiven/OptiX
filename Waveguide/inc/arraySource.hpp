#ifndef ARRAY_SOURCE_H
#define ARRAY_SOURCE_H
#include "fixedValuesSource.hpp"
#include <armadillo>

/** Similar to FixedValuesSource, except that it does not interpolate. Thus, the size of the array must exactly match the simulation grid */
class ArraySource: public FixedValuesSource
{
public:
  ArraySource(): FixedValuesSource("arraySource"){};

  /** Return the value array */
  const arma::cx_vec& getVec() const { return *values; };
};

#endif
