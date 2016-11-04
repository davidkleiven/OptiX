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
  void threePointStencil( unsigned int centerIndx, unsigned int iz, cdouble &left, cdouble &center, cdouble &right ) const;
  double getShiftedX( double x, unsigned int iter ) const;
  void setWG( const WaveGuideFDSimulation &wg );
  std::vector<int> getAccumulatedPixelShift() const { return accumulatedPixelShift; }; // Only for debugging
private:
  const WaveGuideFDSimulation *wg{NULL};
  std::vector<int> accumulatedPixelShift;
  unsigned int border;
  int shiftRelativeToPrevious{0};
};
#endif
