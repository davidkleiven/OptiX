#include "sincSrc.h"
#include <cmath>
#define TINY 1E-6

using namespace std;
const double PI = acos(-1.0);
complex<double> IMAG_UNIT(0.0,1.0);

sinc_src_time::sinc_src_time(double f, double fwidth):freq(f), width(fwidth)
{
  const double zeroAmp = 0.05;
  double domega = 2.0*PI*width;
  peak_time = 2.0/(zeroAmp*width);
};

sinc_src_time::sinc_src_time(const sinc_src_time &copy):freq(copy.freq), width(copy.width), peak_time(copy.peak_time){};

std::complex<double> sinc_src_time::dipole(double time) const
{
  if ( time > 2.0*peak_time )
  {
    return 0.0;
  }

  const double domega = 2.0*PI*width;
  const double omega0 = 2.0*PI*freq;
  double sincArg = domega*(time-peak_time)/2.0;
  double sincVal = 0.0;
  if ( abs(sincArg) < TINY )
  {
    // Use series expansion around zero
    sincVal = 1.0-sincArg*sincArg/6.0;
  }
  else
  {
    sincVal = sin(sincArg)/sincArg;
  } 
  return domega*exp(-IMAG_UNIT*omega0*(time-peak_time))*sincVal;
}

double sinc_src_time::last_time() const
{
  return 2.0*peak_time;
}

bool sinc_src_time::is_equal(const meep::src_time &srcTime) const
{
  const sinc_src_time* other = dynamic_cast<const sinc_src_time*>(&srcTime);
  if ( other )
  {
    return ( this->freq == other->freq ) && ( this->width == other->width ) && ( this->peak_time == other->peak_time );  
  }
  return false;
}

