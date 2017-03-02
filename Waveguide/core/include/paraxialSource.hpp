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
  enum class Dim_t {TWO_D, THREE_D};
  ParaxialSource( const char* name ):name(name){};

  /** Evaluate the source at position x and z */
  virtual cdouble get( double x, double z ) const;

  /** Get function in 3D */
  virtual cdouble get( double x, double y, double z ) const;

  /** Set the wavenumber in nm^{-1} */
  void setWavenumber( double k ){ wavenumber = k; };

  /** Get the wavenumber */
  double getWavenumber() const;

  /** Set the wavelength in nano meters */
  void setWavelength( double lambda );

  /** Set origin of the source */
  void setOrigin( double xOrig, double zOrig ){z0=zOrig;};

  /** Set the overall constant amplitude factor */
  void setAmplitude( double amp ){ amplitude = amp; };

  /** Set the dimension, 2D is default */
  void setDim( Dim_t newdim ){ dim = newdim; };

  /** Get the dimension */
  Dim_t getDim() const { return dim; };

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
  Dim_t dim{Dim_t::TWO_D};
private:
  std::string name;
};
#endif
