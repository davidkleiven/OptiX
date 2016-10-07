#ifndef WAVE_GUIDE_1D_SIMULATION_H
#define WAVE_GUIDE_1D_SIMULATION_H
#include <cstddef>
class Cladding;

class WaveGuide1DSimulation
{
public:
  WaveGuide1DSimulation();
  void setCladding( const Cladding &cladding );
  void setWavenumber( double k ){ wavenumber = k; };
  virtual double potential( double x ) const = 0;
  double getWavenumber() const { return wavenumber; };
protected:
  double width;
  double innerRadius;
  double outerRadius;
  const Cladding *cladding{NULL};
  double wavenumber;
};
#endif
