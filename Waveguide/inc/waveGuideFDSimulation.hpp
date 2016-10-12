#ifndef WAVE_GUIDE_FD_SIMULATION_H
#define WAVE_GUIDE_FD_SIMULATION_H
#include <cstddef>
#include <complex>
#include <string>
#include <jsoncpp/json/writer.h>
class Solver2D;

struct Disctretization
{
  double min;
  double max;
  double step;
};

typedef std::complex<double> cdouble;

class Cladding;
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
  double getWavenumber() const{ return wavenumber; };
  void setWavenumber( double k ){ wavenumber = k; };
  void setWaveLength( double lambda );
  void setSolver( Solver2D &solv );
  std::string getName() const { return name; };
  void setCladding( const Cladding &clad );
  const Cladding& getCladding() const { return *cladding; };
  void solve();
  void save( const std::string &fname ) const;
  void saveWG( const std::string &fname ) const;

  // Virtual methods
  // Refractive index: n = 1 - delta + i*beta
  virtual void getXrayMatProp( double x, double z, double &delta, double &beta ) const = 0;
  virtual void setBoundaryConditions() = 0; // This function should fill the boundary
  virtual void fillInfo( Json::Value &obj ) const {};
  virtual cdouble transverseBC( double z ) const{ return 0.0; };
protected:
  Solver2D *solver{NULL};
  Disctretization *xDisc; // Transverse
  Disctretization *zDisc; // Along optical axis
  double wavenumber;
  std::string name;
  const Cladding *cladding{NULL};

  double* allocateSolutionMatrix() const;
  void deallocateSolutionMatrix( double *matrix ) const;

  void sparseSave( const std::string &fname, double intensityThreshold ) const;
  // Virtual funcitons
  virtual bool isInsideGuide( double x, double z ) const { return true; };
};
#endif
