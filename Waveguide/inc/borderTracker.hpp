#ifndef BORDER_TRACKER_H
#define BORDER_TRACKER_H
#include <complex>
#include <vector>

typedef std::complex<double> cdouble;

class WaveGuideFDSimulation;
class BorderTracker
{
public:
  BorderTracker();
  void init();
  void locateBorder( double z );
  int accumulatedShift() const;
  void threePointStencil( unsigned int centerIndx, double z, cdouble &left, cdouble &center, cdouble &right ) const;
  double getShiftedX( double x ) const;
  void setWG( const WaveGuideFDSimulation &wg );
private:
  const WaveGuideFDSimulation *wg{NULL};
  std::vector<int> accumulatedPixelShift;
  unsigned int border;
  int shiftRelativeToPrevious{0};
  unsigned int currentIteration{0};
};
#endif
