#include "exitFieldSource.hpp"
#include <complex>
#include <stdexcept>

using namespace std;
typedef complex<double> cdouble;

void ExitFieldSource::load( const string &amp, const string &phase )
{
  cdouble im(0.0,1.0);
  arma::vec ph;
  arma::vec ampl;
  ph.load(phase.c_str(), arma::hdf5_binary);
  ampl.load(amp.c_str(), arma::hdf5_binary);
  values = ampl*arma::cos(ph) + im*ampl*arma::sin(ph);
  hasLoadedData = true;
}

void ExitFieldSource::setDiscretization( double min, double max )
{
  if ( !hasLoadedData )
  {
    throw (runtime_error("Need to load data first!"));
  }
  xDisc.min = min;
  xDisc.max = max;
  xDisc.step = (max-min)/values.n_elem;
}

unsigned int ExitFieldSource::closest( double x ) const
{
  return (x-xDisc.min)/xDisc.step;
}

cdouble ExitFieldSource::operator()( double x ) const
{
  if (( x > xDisc.max ) || ( x < xDisc.min ))
  {
    throw (runtime_error("x out of bounds!"));
  }
  unsigned int indx = closest(x);
  double x0 = xDisc.min + indx*xDisc.step;
  double x1 = x0 + xDisc.step;
  return ( (x1-x)*values(indx+1) + (x-x0)*values(indx) )/xDisc.step;
}
