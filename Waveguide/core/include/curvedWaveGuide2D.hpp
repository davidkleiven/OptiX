#ifndef WAVE_GUIDE_CURVED_2D_H
#define WAVE_GUIDE_CURVED_2D_H
#include "waveGuideFDSimulation.hpp"
#include "transmittivity.hpp"
#include <complex>
#include <vector>
#include <armadillo>

class ControlFile;

typedef std::complex<double> cdouble;

/** Class for handling the curved waveguide geometry */
class CurvedWaveGuideFD: public WaveGuideFDSimulation
{
public:
  CurvedWaveGuideFD();
  CurvedWaveGuideFD( const CurvedWaveGuideFD &other );
  CurvedWaveGuideFD operator =(const CurvedWaveGuideFD &other );

  virtual ~CurvedWaveGuideFD();

  /** Set the radius of curvature in nano meters */
  void setRadiusOfCurvature( double newR ) { R = newR; };

  /** Set the width of the wavguide in nano meters */
  void setWidth( double newWidth ) { width = newWidth; };

  /** Get the width of the waveguide in nano meters */
  double getWidth() const { return width; };

  /** Get the radius of curvature in nano meter */
  double getRadiusOfCurvature() const { return R; };

  /** Get the field inside the waveguide */
  void getFieldInsideWG( arma::mat &matrix ) const;

  /** Project the solution onto the eigenmodes */
  double project( double z, const WaveGuide1DSimulation &eig, unsigned int eigenmode ) const;

  /** Returns the weight factor if a smoothed refractive index profile is used for instance to mimic roughness */
  double smoothedWG( double x, double z ) const;

  /** Returns the transmittivity module */
  const post::Transmittivity& getTransmittivity() { return *transmittivity; };

  /** Use a smoothed refractive index profile */
  void useSmoothedWG(){ useSmoothed=true; };

  // Virtual functions
  /** Fill a JSON object with parameters specific to curved waveguides */
  virtual void fillInfo( Json::Value &obj ) const override;

  /** Initialize the simulation. Relevant if the results from an old simulation is loaded to redo the post processing */
  virtual void init( const ControlFile &ctl ) override;

  /** Checks if the coordinates given is inside the waveguide */
  virtual bool isInsideGuide( double x, double z ) const override;

  /** Get material properties. Overloaded in case a smoothed waveguide profile is used. */
  virtual void getXrayMatProp( double x, double z, double &delta, double &beta ) const override;

  /** Extracts the field inside the waveguide at distance wcrd from the edge
  * Thus: wcrd = 0.0: along the lower wall. wcrd=width extracts along the upper wall */
  virtual void extractField( double wcrd, std::vector<cdouble> &res ) const {};
  //virtual void extractField( double wcrd, std::vector<double> &res ) const {} // real part of field

  /** Solve the simulation */
  virtual void solve() override;

  /** Save the results */
  virtual void save( ControlFile &ctl ) override;

  /** Returns the lower border of the waveguide at position z. Required when computing the transmission. */
  virtual double waveGuideStartX( double z ) const;

  /** Returns the upper part of the border at position z. Required when computing the transmission. */
  virtual double waveGuideEndX( double z ) const;
protected:
  CurvedWaveGuideFD( const char *name);
  double R;
  double width;
  bool useSmoothed{false};
  post::Transmittivity *transmittivity{NULL};
};
#endif
