#ifndef SINC_SOURCE_H
#define SINC_SOURCE_H
#include "meep.hpp"
#include <complex>

class sinc_src_time: public meep::src_time
{
public:
  sinc_src_time(double f, double fwidth, double peak_time);
  sinc_src_time(const sinc_src_time &copy);

  virtual std::complex<double> dipole(double t) const;
  virtual double last_time() const;
  virtual src_time* clone() const { return new sinc_src_time(*this);};
  virtual bool is_equal(const meep::src_time &srcTime) const;
  virtual std::complex<double> frequency() const {return freq;};
  virtual void set_frequency(std::complex<double> f){ freq = real(f);};
private:
  double freq;
  double width;
  double peak_time;
};

  
#endif
