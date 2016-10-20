#ifndef WAVE_GUIDE_FD_SIMULATION_H
#define WAVE_GUIDE_FD_SIMULATION_H
#include <cstddef>
#include <complex>
#include <string>
#include <jsoncpp/json/writer.h>
#include <armadillo>
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

class WaveGuideFDSimulation
{
public:
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
  double getWavenumber() const{ return wavenumber; };
  void setWavenumber( double k ){ wavenumber = k; };
  void setWaveLength( double lambda );
  void setSolver( Solver2D &solv );
  std::string getName() const { return name; };
  void setCladding( const Cladding &clad );
  const Cladding& getCladding() const { return *cladding; };
  void solve();
  void save( ControlFile &ctl ) const;
  void save( ControlFile &ctl, double intensityThreshold ) const;
  void saveWG( const std::string &fname ) const;
  double getIntensity( double x, double z ) const; // Using linear interpolation
  double getIntensity( unsigned int ix, unsigned int iz ) const; // Returns value in matrix at (ix,iz)

  // Refractive index: n = 1 - delta + i*beta
  void getXrayMatProp( double x, double z, double &delta, double &beta ) const;

  // Virtual methods
  virtual void setBoundaryConditions() = 0; // This function should fill the boundary
  virtual void fillInfo( Json::Value &obj ) const {};
  virtual cdouble transverseBC( double z ) const{ return 0.0; };
  virtual void init( const ControlFile &ctl );
protected:
  Solver2D *solver{NULL};
  Disctretization *xDisc; // Transverse
  Disctretization *zDisc; // Along optical axis
  double wavenumber;
  std::string name;
  const Cladding *cladding{NULL};
  arma::vec *farFieldModulus{NULL};

  double* allocateSolutionMatrix() const;
  void deallocateSolutionMatrix( double *matrix ) const;
  bool solverInitializedViaInit{false};

  void sparseSave( const std::string &fname, double intensityThreshold ) const;
  void closestIndex( double x, double z, unsigned int &ix, unsigned int &iz ) const;
  double getZ( unsigned int iz ) const;
  double getX ( unsigned int ix ) const;
  double trapezoidalIntegrateIntensityZ( unsigned int iz, unsigned int ixStart, unsigned int ixEnd ) const;
  void getExitField( arma::vec &vec ) const;
  void saveFarField( const std::string &fname, unsigned int uid ) const;

  // Virtual funcitons
  virtual bool isInsideGuide( double x, double z ) const { return true; };
};
#endif
