#ifndef COUPLED_CURVED_WG_H
#define COUPLED_CURVED_WG_H
#include "waveGuideFDSimulation.hpp"

#include <complex>
#include <jsoncpp/json/writer.h>

class CurvedWaveGuideFD;
class ControlFile;

typedef std::complex<double> cdouble;
class CoupledCurvedWG: public WaveGuideFDSimulation
{
public:
  enum class Coordinate_t{CARTESIAN, CYLINDRICAL};
  CoupledCurvedWG(Coordinate_t crdsystem);
  virtual ~CoupledCurvedWG();
  void setSeparation( double sep ) { separation = sep; };
  double getSeparation() const { return separation; };
  void setStartCoupler( double start ){ startCoupler = start; };
  cdouble transverseBC( double z, WaveGuideFDSimulation::Boundary_t bnd) const override;

  CurvedWaveGuideFD& getWg1() { return *wg1; };
  CurvedWaveGuideFD& getWg2() { return *wg2; };

  // Overrides
  void fillInfo( Json::Value &obj ) const override;
  void init( const ControlFile &ctl ) override;
protected:
  Coordinate_t crd;
  double separation;
  double startCoupler;
  CurvedWaveGuideFD *wg1{NULL};
  CurvedWaveGuideFD *wg2{NULL};

  bool isInsideGuide( double x, double z ) const override;
};

#endif
