#include "FluxIntegrator.h"
#include <complex>
using namespace std;

FluxIntegrator::FluxIntegrator(): _zpos(0.0), _nEvalPointsInEachDirection(0), _evaluationCrd(NULL), _flux(NULL),\
_xmin(0.0), _xmax(1.0), _ymin(0.0), _ymax(1.0){};

FluxIntegrator::~FluxIntegrator()
{
  if ( _evaluationCrd != NULL )
  {
    delete _evaluationCrd;
  }

  if ( _flux != NULL )
  {
    delete _flux;
  }
}

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

void FluxIntegrator::fillEvaluationXY()
{
  if ( _nEvalPointsInEachDirection == 0 )
  {
    this->setnEvalPoints(20);
  } 

  // Fill evaluation points
  double dx = (_xmax-_xmin)/static_cast<double>(_nEvalPointsInEachDirection);
  double dy = (_ymax-_ymin)/static_cast<double>(_nEvalPointsInEachDirection);
  double x = _xmin;
  for ( unsigned int ix=0;ix<_nEvalPointsInEachDirection;ix++ )
  {
    double y = _ymin;
    for ( unsigned int iy=0;iy<_nEvalPointsInEachDirection;iy++)
    {
      _evaluationCrd->SetEntry( ix*_nEvalPointsInEachDirection+iy, 0, x);
      _evaluationCrd->SetEntry( ix*_nEvalPointsInEachDirection+iy, 1, y);
      _evaluationCrd->SetEntry( ix*_nEvalPointsInEachDirection+iy, 2, _zpos);
      y += dy;
    }
    x += dx;
  }
}    

double FluxIntegrator::compute( scuff::RWGGeometry &geo, IncField *IF, HVector *vec, double omega, double kBloch[2] )
{
  if ( _nEvalPointsInEachDirection == 0 )
  {
    this->fillEvaluationXY();
  }

  // Get field at all evaluation points      
  geo.GetFields(IF, vec, omega, kBloch, _evaluationCrd, _flux);
  complex<double> Ex, Ey, Hx, Hy; 
  double zFlux = 0.0;
  for ( unsigned int i=0;i<_nEvalPointsInEachDirection;i++ )
  {
    Ex = _flux->GetEntry(i,0);
    Ey = _flux->GetEntry(i,1);
    Hx = _flux->GetEntry(i,3);
    Hy = _flux->GetEntry(i,4);
    zFlux += 0.5*real( Ex*conj(Hy) - conj(Ey)*Hx );
  }
  return zFlux/static_cast<double>( _nEvalPointsInEachDirection*_nEvalPointsInEachDirection );
}
    
double FluxIntegrator::incidentFlux( scuff::RWGGeometry &geo, IncField &IF, double omega, double kBloch[2] )
{
  return this->compute( geo, &IF, NULL, omega, kBloch );
}

double FluxIntegrator::scatteredFlux( scuff::RWGGeometry &geo, HVector &vec, double omega, double kBloch[2] )
{
  return this->compute( geo, NULL, &vec, omega, kBloch );
}
