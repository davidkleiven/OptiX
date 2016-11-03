#include "borderTracker.hpp"
#include "WaveGuideFDSimulation.hpp"
#include <stdexcept>

using namespace std;

BorderTracker::BorderTracker()
{
  accumulatedShift.push_back(0);
}

void BorderTracker::locateBorder( double z )
{
  if ( wg == NULL )
  {
    throw ( runtime_error( "BorderTracker: No waveguide is set!") );
  }

  bool borderFound = false;
  unsigned int shift = 0;
  int prevShift = accumulatedShift[ accumulatedShift.size()-1 ];
  while ( !borderFound )
  {
    double xBorder = wg->getX( upperBorder+shift );
    double xNeighbour = wg->getX( upperBorder+shift+1 );
    if ( isInside(xBorder,z) && !isInside(xNeighbour,z) )
    {
      upperBorder += shift;
      pixelShift.push_back( prevShift + shift);
      borderFound = true;
    }

    // Try the other direction
    if ( upperBorder < shift ) continue;


    xBorder = wg->getX( upperBorder - shift );
    xNeighbour = wg->getX( upperBorder - shift + 1);
    if ( isInside(xBorder, z) && !isInside(xNeighbour,z) )
    {
      upperBorder -= shift;
      pixelShift.push_back( prevShift-static_cast<int>(shift) );
      borderFound = true;
    }
    shift++;

    if ( upperBorder + shift  > wg->getSolver().getSolution().n_rows )
    {
      throw ( runtime_error("BorderTracker: Could not locate border!") );
    }
  }
}

double BorderTracker::getShiftedX( double x ) const
{
  double stepX = wg->transverseDiscretization().step;
  return x + accumulatedShift[accumulatedShift.size()-1]*stepX;
}
