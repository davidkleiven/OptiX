#ifndef COUPLED_CURVED_WG_H
#define COUPLED_CURVED_WG_H
#include "waveGuideFDSimulation.hpp"

#include <complex>
#include <jsoncpp/json/writer.h>

class CurvedWaveGuideFD;
class ControlFile;

typedef std::complex<double> cdouble;

/** Class for handling the Coupled curved waveguide geometry */
class CoupledCurvedWG: public WaveGuideFDSimulation
{
public:
  enum class Coordinate_t{CARTESIAN, CYLINDRICAL};
  CoupledCurvedWG(Coordinate_t crdsystem);
  virtual ~CoupledCurvedWG();

  /** Set the radial separation between the waveguides in nano meters */
  void setSeparation( double sep ) { separation = sep; };

  /** Get the separation between the waveguides in nano meters */
  double getSeparation() const { return separation; };

  /** Set the position where the waveguide coupling to the other starts. This can be used to ensure that the direct beam only enters one of the waveguides */
  void setStartCoupler( double start ){ startCoupler = start; };

  /** Specifies the boundary conditions at x = +-L, where the simulation domain is -L < x < L */
  cdouble transverseBC( double z, WaveGuideFDSimulation::Boundary_t bnd) const override;

  /** Get the first waveguide */
  CurvedWaveGuideFD& getWg1() { return *wg1; };

  /** Get the second waveguide (the coupler) */
  CurvedWaveGuideFD& getWg2() { return *wg2; };

  // Overrides
  /** Fill a JSON object with parameters specific to this geometry */
  void fillInfo( Json::Value &obj ) const override;

  /** Initialize the geometry. Only when loading an old simulation for the purpose of running the postprocessing */
  void init( const ControlFile &ctl ) override;
protected:
  Coordinate_t crd;
  double separation;
  double startCoupler;
  CurvedWaveGuideFD *wg1{NULL};
  CurvedWaveGuideFD *wg2{NULL};

  /** Checks if the given coordinates are inside the waveguide */
  bool isInsideGuide( double x, double z ) const override;
};

#endif
