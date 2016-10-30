#ifndef EXIT_FIELD_SOURCE_H
#define EXIT_FIELD_SOURCE_H
#include <armadillo>
#include <string>
#include "waveGuideFDSimulation.hpp"

class ExitFieldSource
{
public:
  ExitFieldSource(){};
  cdouble operator()(double x) const;
  void setDiscretization( double xmin, double xmax );
  void load( const std::string &amp, const std::string &phase );
private:
  arma::cx_vec values;
  Disctretization xDisc;
  unsigned int closest( double x ) const;
  bool hasLoadedData{false};
};
#endif
