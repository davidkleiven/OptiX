#ifndef BEAM_SPLITTER_H
#define BEAM_SPLITTER_H
#include "waveGuideFDSimulation.hpp"
#include <jsoncpp/json/writer.h>

class BeamSplitter: public WaveGuideFDSimulation
{
public:
  BeamSplitter(): WaveGuideFDSimulation("BeamSplitter"){};
  void setAngleWithZAxisDeg( double angle ){ angleWithZAxisDeg=angle; };
  void setSplitStart( double start ){splitStart=start;};
  void setWidth( double w ){ width=w; };

  // Overrides
  virtual void fillInfo( Json::Value &obj ) const override final;
private:
  double angleWithZAxisDeg{0.0};
  double splitStart{0.0};
  double width{100.0};

  //Overrides
  bool isInsideGuide( double x, double z ) const override;
};
#endif
