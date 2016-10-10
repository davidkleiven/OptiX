#ifndef WAVE_GUIDE_FD_SIMULATION_H
#define WAVE_GUIDE_FD_SIMULATION_H
#include <cstddef>
#include <complex>
class Solver2D;

struct Disctretization
{
  double min;
  double max;
  double step;
};

class WaveGuideFDSimulation
{
public:
  WaveGuideFDSimulation();
  ~WaveGuideFDSimulation();
  void setTransverseDiscretization( double xmin , double xmax, double step );
  void setLongitudinalDiscretization( double zmin , double zmax, double step );
  unsigned int nodeNumberTransverse() const;
  unsigned int nodeNumberLongitudinal() const;
  const Disctretization& transverseDiscretization() const{ return *xDisc; };
  const Disctretization& longitudinalDiscretization() const { return *zDisc; };
  virtual std::complex<double> getRefractiveIndex( double x, double z ) const = 0;
  double getWavenumber() const{ return wavenumber; };
  void setWavenumber( double k ){ wavenumber = k; };
private:
  Solver2D *solver{NULL};
  Disctretization *xDisc; // Transverse
  Disctretization *zDisc; // Along optical axis
  double wavenumber;
};
#endif
