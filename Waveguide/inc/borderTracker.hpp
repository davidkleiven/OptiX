#ifndef BORDER_TRACKER_H
#define BORDER_TRACKER_H

class WaveGuideFDSimulation;
class BorderTracker
{
public:
  BorderTracker();
  void locateBorder( double z );
  int accumulatedShift() const;
  double getShiftedX( double x ) const;
private:
  const WaveGuideFDSimulation *wg{NULL};
  std::vector<int> accumulatedPixelShift;
  unsigned int upperBorder;
};
#endif
