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

  /** Evaluates the field at position x */
  cdouble operator()(double x) const;

  /** Sets the transverse discretization */
  void setDiscretization( double xmin, double xmax );

  /** Returns the transverse discretization */
  Disctretization getDisc() const { return xDisc; };

  /** Load field from HDF5 file */
  void load( H5::H5File &file );

  /** Not used */
  cdouble initfield( double x ) const;
private:
  arma::cx_vec values;
  Disctretization xDisc;
  unsigned int closest( double x ) const;
  bool hasLoadedData{false};
};
#endif
