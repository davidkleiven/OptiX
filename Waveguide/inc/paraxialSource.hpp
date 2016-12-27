#ifndef SOURCE_H
#define SOURCE_H
#include <complex>
#include <string>
#include <json/writer.h>

typedef std::complex<double> cdouble;

/** Base class for the paraxial sources */
class ParaxialSource
{
public:
  ParaxialSource( const char* name ):name(name){};

  /** Evaluate the source at position x and z */
  virtual cdouble get( double x, double z ) const = 0;

  /** Set the wavenumber in nm^{-1} */
  void setWavenumber( double k ){ wavenumber = k; };

  /** Set the wavelength in nano meters */
  void setWavelength( double lambda );

  /** Set origin of the source */
  void setOrigin( double xOrig, double zOrig ){z0=zOrig;};

  /** Set the overall constant amplitude factor */
  void setAmplitude( double amp ){ amplitude = amp; };

  /** Get the name of the source */
  std::string getName() const { return name; };

  /** Get the wavelength in nano meters */
  double getWavelength() const;

  /** Fill JSON object with parameters specific to this class */
  virtual void info( Json::Value &obj ) const;
protected:
  double wavenumber{0.0};
  double z0{0.0};
  double x0{0.0};
  double amplitude{1.0};
private:
  std::string name;
};
#endif
