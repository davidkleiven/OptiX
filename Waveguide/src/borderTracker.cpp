#include "borderTracker.hpp"
#include "waveGuideFDSimulation.hpp"
#include "solver2D.hpp"
#include <stdexcept>
#include <cassert>

using namespace std;

BorderTracker::BorderTracker()
{
  accumulatedPixelShift.push_back(0);
}

void BorderTracker::locateBorder( double z )
{
  assert ( wg != NULL );

  bool borderFound = false;
  unsigned int shift = 0;
  int prevShift = accumulatedPixelShift[ accumulatedPixelShift.size()-1 ];
  while ( !borderFound )
  {
    double xBorder = wg->getX( border+shift );
    double xNeighbour = wg->getX( border+shift+1 );
    if ( wg->isInsideGuide(xBorder,z) && !wg->isInsideGuide(xNeighbour,z) )
    {
      border += shift;
      accumulatedPixelShift.push_back( prevShift + shift);
      shiftRelativeToPrevious = shift;
      borderFound = true;
    }

    // Try the other direction
    if ( border < shift ) continue;


    xBorder = wg->getX( border - shift );
    xNeighbour = wg->getX( border - shift + 1);
    if ( wg->isInsideGuide(xBorder, z) && !wg->isInsideGuide(xNeighbour,z) )
    {
      border -= shift;
      accumulatedPixelShift.push_back( prevShift-static_cast<int>(shift) );
      shiftRelativeToPrevious = -static_cast<int>(shift);
      borderFound = true;
    }
    shift++;

    if (( border + shift  >= wg->getSolver().getSolution().n_rows ) || ( border-shift < 0 ))
    {
      throw ( runtime_error("BorderTracker: Could not locate border!") );
    }
  }
  currentIteration++;
}

double BorderTracker::getShiftedX( double x ) const
{
  assert ( wg != NULL );
  double stepX = wg->transverseDiscretization().step;
  return x + accumulatedPixelShift[accumulatedPixelShift.size()-1]*stepX;
}

void BorderTracker::threePointStencil( unsigned int centerIndx, double z, cdouble &left, cdouble &center, cdouble &right ) const
{
  assert( wg != NULL );
  int indx = centerIndx + shiftRelativeToPrevious;
  unsigned int maxIndx = wg->getSolver().getSolution().n_rows;

  if ( ( indx-1 < 0 ) || ( indx-1 > maxIndx ) )
  {
    left = wg->transverseBC(z);
  }
  else
  {
    left = wg->getSolver().getSolution()( indx-1, currentIteration-1);
  }

  if (( indx < 0 ) || ( indx > maxIndx ))
  {
    center = wg->transverseBC(z);
  }
  else
  {
    center = wg->getSolver().getSolution()( indx, currentIteration-1);
  }

  if (( indx +1 < 0 ) || ( indx+1 > maxIndx ))
  {
    right = wg->transverseBC(z);
  }
  else
  {
    right = wg->getSolver().getSolution()(indx+1, currentIteration-1);
  }
}

void BorderTracker::setWG( const WaveGuideFDSimulation &newGuide )
{
  wg = &newGuide;
}

void BorderTracker::init()
{
  assert( wg != NULL );

  // TODO: Should be generalized to z != 0. Need to be done when generalizing the source to be at z != 0
  double z = 0.0;
  bool borderFound = false;
  // Find the border
  for ( unsigned int ix=0;ix<wg->nodeNumberTransverse()-1;ix++ )
  {
    double x = wg->getX(ix);
    double xNext = wg->getX(ix+1);
    if ( wg->isInsideGuide(x,z) && !wg->isInsideGuide(xNext,z) )
    {
      borderFound = true;
      break;
    }
  }

  if ( !borderFound )
  {
    throw ( runtime_error("BorderTracker: Did no manage to locate the border during initialization!") );
  }
}
