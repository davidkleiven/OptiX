#ifndef WAVE_GUIDE_1D_SIMULATION_H
#define WAVE_GUIDE_1D_SIMULATION_H
#include <cstddef>
#include <string>
#include <complex>

typedef std::complex<double> cdouble;
class Cladding;
class Solver1D;
class ControlFile;

/** Base class for 1D waveguides */
class WaveGuide1DSimulation
{
public:
  WaveGuide1DSimulation(const char* name):name(name){};
  virtual ~WaveGuide1DSimulation();

  /** Set cladding */
  void setCladding( const Cladding &cladding );

  /** Set material inside */
  void setCore( const Cladding &core );

  /** Get cladding */
  const Cladding& getCladding() const { return *cladding; };

  /** Set waveguide width in nano meters */
  void setWidth( double newwidth ){ width = newwidth; };

  /** Get waveguide width in nano meters */
  double getWidth() const { return width; };

  /** Set 1D solver */
  void setSolver( Solver1D &solv );

  /** Set wavenumber in nm^{-1} */
  void setWavenumber( double k ){ wavenumber = k; };

  /** Set wavelength in nano meters */
  void setWaveLength( double lambda ){ wavenumber = 2.0*PI/lambda; };

  /** Get induced potential by the refractive index in nm^{-2} */
  virtual double potential( double x ) const {return 0.0;};

  /** Get complex potential in nm^{-2} */
  virtual cdouble complexPotential( double x ) const { return 0.0; };

  /** Get wavenumber in nm^{-1} */
  double getWavenumber() const { return wavenumber; };

  /** Run the simulation */
  void solve();

  /** Save results to HDF5 files */
  void save( ControlFile &ctl ) const;

  /** Load an old solution */
  virtual void load( ControlFile &ctl );

  /** Get solver */
  const Solver1D* getSolver() const { return solver; };

  /** Use complex potential */
  bool useComplex() const { return useComplexPotential; };
protected:
  double width;
  double innerRadius;
  double outerRadius;
  const Cladding *cladding{NULL};
  const Cladding *core{NULL};
  double wavenumber;
  Solver1D *solver{NULL};
  static const double PI;
  std::string name;

  /** Save potential to file */
  void writePotentialToFile( const std::string &fname, double xmin, double xmax ) const;
  bool useComplexPotential{false};
private:
    bool solverInitializedFromLoad{false};
};
#endif
