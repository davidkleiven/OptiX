#ifndef HARMONIC_INVERSION_WRAPPER_H
#define HARMONIC_INVERSION_WRAPPER_H
#include <vector>
#include <complex>
#include <map>
#include <string>

typedef std::complex<double> cdouble;

struct HarmInvRes
{
  std::vector<double> freq;
  std::vector<double> amplitude;
  std::vector<double> phase;
  std::vector<double> data;
  std::vector<double> decay;
  std::map<std::string, double> attrib;
};

class HarmonicInversion
{
public:
  HarmonicInversion(){};
  void solve( std::vector<cdouble> &data );
  void setFreq( double freqMin, double freqMax, unsigned int nFreq );
  void addAttribute( const char* name, double value, unsigned int dset );
  void save( const char* fname ) const;
private:
  std::vector<HarmInvRes> results;
  double freqMin{0.0};
  double freqMax{0.0};
  unsigned int nFreq{0};
};

#endif
