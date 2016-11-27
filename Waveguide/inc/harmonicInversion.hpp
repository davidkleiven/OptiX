#ifndef HARMONIC_INVERSION_WRAPPER_H
#define HARMONIC_INVERSION_WRAPPER_H
#include <vector>
#include <complex>
#include <map>
#include <string>

typedef std::complex<double> cdouble;

/** Struct for storing the output result from the Harminv library*/
struct HarmInvRes
{
  std::vector<double> freq;
  std::vector<double> amplitude;
  std::vector<double> phase;
  std::vector<double> data;
  std::vector<double> decay;
  std::map<std::string, double> attrib;
};


/** Wrapper class for the Harminv library */
class HarmonicInversion
{
public:
  HarmonicInversion(){};
  /** Compute the harmonic inversion */
  void solve( std::vector<cdouble> &data );

  /** Set the frequency domain */
  void setFreq( double freqMin, double freqMax, unsigned int nFreq );

  /** Add attribute to dataset number given by dset */
  void addAttribute( const char* name, double value, unsigned int dset );

  /** Save the results to a HDF5 file */
  void save( const char* fname ) const;
private:
  std::vector<HarmInvRes> results;
  double freqMin{0.0};
  double freqMax{0.0};
  unsigned int nFreq{0};
};

#endif
