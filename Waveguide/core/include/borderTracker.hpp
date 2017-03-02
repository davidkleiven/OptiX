#ifndef BORDER_TRACKER_H
#define BORDER_TRACKER_H
#include <complex>
#include <vector>

typedef std::complex<double> cdouble;

class WaveGuideFDSimulation;

/** Class that tracks the border of a geometry */
class BorderTracker
{
public:
  BorderTracker();

  /** Initialize the tracker. This function has to be called before anything else */
  void init();

  /** Locates the border at position z*/
  void locateBorder( double z );

  /** Returns the total number of indices the border has shifted */
  int accumulatedShift() const;

  /** Returns the three point stencil required in the 1D Cranck-Nicholson FD-scheme */
  void threePointStencil( unsigned int centerIndx, unsigned int iz, cdouble &left, cdouble &center, cdouble &right ) const;

  /** Get the x coordinate corresponding to the old after the shift */
  double getShiftedX( double x, unsigned int iter ) const;

  /** Set the waveguide for the BorderTracker to operate on */
  void setWG( const WaveGuideFDSimulation &wg );

  /** Get the accumaleted shift for each index of the matrix */
  std::vector<int> getAccumulatedPixelShift() const { return accumulatedPixelShift; }; // Only for debugging
private:
  const WaveGuideFDSimulation *wg{NULL};
  std::vector<int> accumulatedPixelShift;
  unsigned int border;
  int shiftRelativeToPrevious{0};
};
#endif
