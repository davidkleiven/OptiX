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
  CoupledCurvedWG();
  virtual ~CoupledCurvedWG();
  void setSeparation( double sep ) { separation = sep; };
  double getSeparation() const { return separation; };
  void setStartCoupler( double start ){ startCoupler = start; };

  CurvedWaveGuideFD& getWg1() { return *wg1; };
  CurvedWaveGuideFD& getWg2() { return *wg2; };

  // Overrides
  void fillInfo( Json::Value &obj ) const override;
  void init( const ControlFile &ctl ) override;
protected:
  double separation;
  double startCoupler;
  CurvedWaveGuideFD *wg1;
  CurvedWaveGuideFD *wg2;

  bool isInsideGuide( double x, double z ) const override;
};

#endif
