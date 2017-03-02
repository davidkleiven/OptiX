#include "borderTracker.hpp"
#include "waveGuideFDSimulation.hpp"
#include "solver2D.hpp"
#include <stdexcept>
#include <cassert>
#include <iostream>

#define DEBUG_BORDER_LOCATE

using namespace std;

BorderTracker::BorderTracker()
{
  accumulatedPixelShift.push_back(0);
}

void BorderTracker::locateBorder( double z )
{
  assert ( wg != NULL );

  bool borderFound = false;
  int shift = 0;
  int prevShift = accumulatedPixelShift.back();
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
      break;
    }

    // Try the other direction
    xBorder = wg->getX( border - shift );
    xNeighbour = wg->getX( border - shift + 1);
    if ( wg->isInsideGuide(xBorder, z) && !wg->isInsideGuide(xNeighbour,z) )
    {
      border -= shift;
      accumulatedPixelShift.push_back( prevShift-static_cast<int>(shift) );
      shiftRelativeToPrevious = -static_cast<int>(shift);
      borderFound = true;
      break;
    }
    shift++;


    if ( shift > wg->getSolver().getSolution().n_rows )
    {
      throw( runtime_error("BorderTracker: The shift is too large!") );
    }
  }
}

double BorderTracker::getShiftedX( double x, unsigned int iter ) const
{
  iter = iter > accumulatedPixelShift.size() ? accumulatedPixelShift.size()-1:iter;
  double stepX = wg->transverseDiscretization().step;
  return x + accumulatedPixelShift[iter]*stepX;
}

void BorderTracker::threePointStencil( unsigned int centerIndx, unsigned int iz, cdouble &left, cdouble &center, cdouble &right ) const
{
  assert( wg != NULL );
  int indx = static_cast<int>(centerIndx) + shiftRelativeToPrevious;
  double z = wg->getZ(iz);
  unsigned int maxIndx = wg->getSolver().getSolution().n_rows-1;
  //cerr << currentIteration << " " << indx << " " << accumulatedPixelShift.back() << endl;
  if ( ( indx-1 < 0 ) || ( indx-1 > maxIndx ) )
  {
    left = 0.0;
  }
  else
  {
    left = wg->getSolver().getSolution()( indx-1, iz-1);
  }
  //cerr << "left\n";
  if (( indx < 0 ) || ( indx > maxIndx ))
  {
    center = 0.0;
  }
  else
  {
    center = wg->getSolver().getSolution()( indx, iz-1);
  }
  //cerr << "center\n";
  if (( indx +1 < 0 ) || ( indx+1 > maxIndx ))
  {
    right = 0.0;
  }
  else
  {
    right = wg->getSolver().getSolution()(indx+1, iz-1);
  }
  //cerr << "right\n";
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
      border = ix;
      borderFound = true;
      break;
    }
  }
  cout << border << endl;
  if ( !borderFound )
  {
    throw ( runtime_error("BorderTracker: Did not manage to locate the border during initialization!") );
  }
}
