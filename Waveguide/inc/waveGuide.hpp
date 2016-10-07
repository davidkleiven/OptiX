#ifndef WAVE_GUIDE_1D_SIMULATION_H
#define WAVE_GUIDE_1D_SIMULATION_H
#include <cstddef>
#include <string>

class Cladding;
class Solver1D;

class WaveGuide1DSimulation
{
public:
  WaveGuide1DSimulation(){};
  void setCladding( const Cladding &cladding );
  void setWidth( double newwidth ){ width = newwidth; };
  void setSolver( Solver1D &solv );
  void setWavenumber( double k ){ wavenumber = k; };
  void setWaveLength( double lambda ){ wavenumber = 2.0*PI/lambda; };
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
  static const double PI;
};
#endif
