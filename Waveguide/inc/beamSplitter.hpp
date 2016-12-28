#ifndef BEAM_SPLITTER_H
#define BEAM_SPLITTER_H
#include "waveGuideFDSimulation.hpp"
#include <json/writer.h>

/** Class handling the beam splitter geometry */
class BeamSplitter: public WaveGuideFDSimulation
{
public:
  BeamSplitter(): WaveGuideFDSimulation("BeamSplitter"){};

  /** Set the angle each of the arms forms with the z-axis. Thus, the angle between the arms is 2 times this angle*/
  void setAngleWithZAxisDeg( double angle ){ angleWithZAxisDeg=angle; };

  /** Set the z-coordinate where the split starts*/
  void setSplitStart( double start ){splitStart=start;};

  /** Set the width of the waveguide in nano meters*/
  void setWidth( double w ){ width=w; };

  // Overrides
  /** Fills a json object with simulation information specific to the BeamSplitter class*/
  virtual void fillInfo( Json::Value &obj ) const override final;
private:
  double angleWithZAxisDeg{0.0};
  double splitStart{0.0};
  double width{100.0};

  //Overrides
  /** Checks if the given coordinates is inside the waveguide */
  bool isInsideGuide( double x, double z ) const override;
};
#endif
