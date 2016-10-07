#ifndef WAVE_GUIDE_1D_SIMULATION_H
#define WAVE_GUIDE_1D_SIMULATION_H
#include <cstddef>
#include <string>
class Cladding;
class Solver1D;

class WaveGuide1DSimulation
{
public:
  WaveGuide1DSimulation();
  void setCladding( const Cladding &cladding );
  void setSolver( Solver1D &solv );
  void setWavenumber( double k ){ wavenumber = k; };
  virtual double potential( double x ) const = 0;
  double getWavenumber() const { return wavenumber; };
  void solve();
  void save( const std::string &fname ) const;
protected:
  double width;
  double innerRadius;
  double outerRadius;
  const Cladding *cladding{NULL};
  double wavenumber;
  Solver1D *solver{NULL};
};
#endif
