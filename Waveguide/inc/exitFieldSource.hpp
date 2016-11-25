#ifndef EXIT_FIELD_SOURCE_H
#define EXIT_FIELD_SOURCE_H
#include <armadillo>
#include <string>
#include "waveGuideFDSimulation.hpp"
#include <H5Cpp.h>

class ExitFieldSource
{
public:
  ExitFieldSource(){};
  cdouble operator()(double x) const;
  void setDiscretization( double xmin, double xmax );
  Disctretization getDisc() const { return xDisc; };
  void load( H5::H5File &file );
  cdouble initfield( double x ) const;
private:
  arma::cx_vec values;
  Disctretization xDisc;
  unsigned int closest( double x ) const;
  bool hasLoadedData{false};
};
#endif
