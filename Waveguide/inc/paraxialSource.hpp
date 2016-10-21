#ifndef SOURCE_H
#define SOURCE_H
#include <complex>
#include <string>

typedef std::complex<double> cdouble;
class ParaxialSource
{
public:
  ParaxialSource( const char* name ):name(name){};
  virtual cdouble get( double x, double z ) const = 0;
  void setWavenumber( double k ){ wavenumber = k; };
  void setWavelength( double lambda );
  void setOrigin( double zOrig ){z0=zOrig;};
  void setAmplitude( double amp ){ amplitude = amp; };
  std::string getName() const { return name; };

  double getWavelength() const;
protected:
    double wavenumber{0.0};
    double z0{0.0};
    double amplitude{1.0};
private:
  std::string name;
};
#endif
