#include "FluxIntegrator.h"

FluxIntegrator::FluxIntegrator(): _zpos(0.0), _nEvalPointsInEachDirection(0), _evaluationCrd(NULL), _flux(NULL){};

void FluxIntegrator::setnEvalPoints(unsigned int nPoints)
{
  if ( _evaluationCrd != NULL )
  {
    delete _evaluationCrd;
  }

  if ( _flux != NULL )
  {
    delete _flux;
  }

  _nEvalPointsInEachDirection = nPoints;
  
  _evaluationCrd = new HMatrix( nPoints*nPoints, 3 );
  _flux = new HMatrix ( nPoints*nPoints, 3 );
}


