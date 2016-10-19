#ifndef WAVE_GUIDE_1D_SIMULATION_H
#define WAVE_GUIDE_1D_SIMULATION_H
#include <cstddef>
#include <string>

class Cladding;
class Solver1D;
class ControlFile;

class WaveGuide1DSimulation
{
public:
  WaveGuide1DSimulation(const char* name):name(name){};
  virtual ~WaveGuide1DSimulation();
  void setCladding( const Cladding &cladding );
  const Cladding& getCladding() const { return *cladding; };
  void setWidth( double newwidth ){ width = newwidth; };
  void setSolver( Solver1D &solv );
  void setWavenumber( double k ){ wavenumber = k; };
  void setWaveLength( double lambda ){ wavenumber = 2.0*PI/lambda; };
  virtual double potential( double x ) const = 0;
  double getWavenumber() const { return wavenumber; };
  void solve();
  void save( ControlFile &ctl ) const;
  void load( ControlFile &ctl );
protected:
  double width;
  double innerRadius;
  double outerRadius;
  const Cladding *cladding{NULL};
  double wavenumber;
  Solver1D *solver{NULL};
  static const double PI;
  std::string name;
  void writePotentialToFile( const std::string &fname, double xmin, double xmax ) const;

private:
    bool solverInitializedFromLoad{false};
};
#endif
