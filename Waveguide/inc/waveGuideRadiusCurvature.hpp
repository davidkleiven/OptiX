#ifndef WAVE_GUIDE_HUGE_CURVATURE_H
#define WAVE_GUIDE_HUGE_CURVATURE_H
#include "waveGuide.hpp"
#include "controlFile.hpp"

/** Class for simulating a waveguide with large radius of curvature relative to its width */
class WaveGuideLargeCurvature: public WaveGuide1DSimulation
{
public:
  WaveGuideLargeCurvature():WaveGuide1DSimulation("LargeRadiusOfCurvature"){};

  /** Set the radius of curvature */
  void setRadiusOfCurvature( double R ){ outerRadius = R; };

  /** Get the potential induced by the refractive index profile */
  double potential( double x ) const override;

  /** Load old solution */
  void load( ControlFile &ctl ) override;
protected:
  WaveGuideLargeCurvature( const char* name ):WaveGuide1DSimulation(name){};
};
#endif
