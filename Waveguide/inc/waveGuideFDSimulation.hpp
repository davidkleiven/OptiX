#ifndef WAVE_GUIDE_FD_SIMULATION_H
#define WAVE_GUIDE_FD_SIMULATION_H
#include <cstddef>
#include <complex>
#include <string>
#include <jsoncpp/json/writer.h>
#include <armadillo>
#include "waveGuideBorder.hpp"
#include "borderTracker.hpp"
class Solver2D;

struct Disctretization
{
  double min;
  double max;
  double step;
};

typedef std::complex<double> cdouble;

class Cladding;
class ControlFile;
class WaveGuide1DSimulation;
class ParaxialSource;

class WaveGuideFDSimulation
{
public:
  enum class Boundary_t {TOP, BOTTOM, LEFT, RIGHT};
  WaveGuideFDSimulation();
  WaveGuideFDSimulation(const char *name);
  ~WaveGuideFDSimulation();
  void setTransverseDiscretization( double xmin , double xmax, double step );
  void setLongitudinalDiscretization( double zmin , double zmax, double step );
  unsigned int nodeNumberTransverse() const;
  unsigned int nodeNumberLongitudinal() const;
  const Disctretization& transverseDiscretization() const{ return *xDisc; };
  const Disctretization& longitudinalDiscretization() const { return *zDisc; };
  void computeFarField();
  void computeFarField( unsigned int signalLength ); // With padding to increase low freq resolution
  double getWavenumber() const{ return wavenumber; };
  void setWavenumber( double k ){ wavenumber = k; };
  void setWaveLength( double lambda );
  void setSolver( Solver2D &solv );
  void setWaveguideLength( double wgl ){ wglength = wgl; };
  std::string getName() const { return name; };
  void setCladding( const Cladding &clad );
  void setInsideMaterial( const Cladding &mat );
  const Cladding& getCladding() const { return *cladding; };
  void solve();
  void save( ControlFile &ctl ) const;
  void save( ControlFile &ctl, double intensityThreshold ) const;
  void saveWG( const std::string &fname ) const;
  void extractWGBorders();
  double getIntensity( double x, double z ) const; // Using linear interpolation
  double getIntensity( unsigned int ix, unsigned int iz ) const; // Returns value in matrix at (ix,iz)
  double getZ( unsigned int iz ) const;
  double getX ( int ix ) const;
  const arma::vec& getFarField() const { return *farFieldModulus; };
  void getExitField( arma::vec &vec ) const;
  const Solver2D& getSolver() const { return *solver; };
  void useBorderTracker();
  BorderTracker* getBorderTracker(){ return bTracker; };
  void closestIndex( double x, double z, unsigned int &ix, unsigned int &iz ) const;

  // Virtual methods
  virtual void setBoundaryConditions( const ParaxialSource& src ); // This function should fill the boundary
  virtual void fillInfo( Json::Value &obj ) const {};
  virtual void init( const ControlFile &ctl );
  virtual cdouble transverseBC( double z, Boundary_t bnd ) const;
  virtual cdouble transverseBC( double z ) const;
  virtual bool isInsideGuide( double x, double z ) const { return true; };
  // Refractive index: n = 1 - delta + i*beta
  virtual void getXrayMatProp( double x, double z, double &delta, double &beta ) const;
protected:
  Solver2D *solver{NULL};
  Disctretization *xDisc; // Transverse
  Disctretization *zDisc; // Along optical axis
  double wavenumber;
  double wglength{1E10};
  std::string name;
  const Cladding *cladding{NULL};
  const Cladding *insideMaterial{NULL};
  arma::vec *farFieldModulus{NULL};
  BorderTracker *bTracker{NULL};

  double* allocateSolutionMatrix() const;
  void deallocateSolutionMatrix( double *matrix ) const;
  bool solverInitializedViaInit{false};

  void sparseSave( const std::string &fname, double intensityThreshold ) const;
  double trapezoidalIntegrateIntensityZ( unsigned int iz, unsigned int ixStart, unsigned int ixEnd ) const;

  void getExitField( arma::cx_vec &vec ) const;
  void saveFarField( const std::string &fname, unsigned int uid ) const;
  const ParaxialSource* src{NULL};
  std::vector<WaveGuideBorder> *wgborder{NULL};

  // Virtual funcitons
  virtual bool waveguideEnded( double x, double z ) const { return z > wglength; };
};
#endif
