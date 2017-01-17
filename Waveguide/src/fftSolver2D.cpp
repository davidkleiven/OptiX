#include "fftSolver2D.hpp"
#include <cmath>
#include "paraxialSimulation.hpp"

const double PI = acos(-1.0);

FFTSolver2D::~FFTSolver2D()
{
  if ( filterSignal != NULL ) delete filterSignal;
}

void FFTSolver2D::init()
{
  initialLength = prevSolution->n_elem;
  unsigned int powerOf2 = log2( initialLength );
  unsigned int newLength = pow( 2, powerOf2+2 );

  arma::cx_vec *ptr = new arma::cx_vec( newLength, 0.0 );
  ptr->subvec( 0, initialLength ) = *prevSolution;
  delete prevSolution;
  prevSolution = ptr;

  delete currentSolution;
  currentSolution = new arma::cx_vec( newLength );
}

cdouble FFTSolver2D::kernel( double kx ) const
{
  cdouble A( 0.0, 0.5/wavenumber );
  return exp( -A*kx*kx );
}

void FFTSolver2D::propagate()
{
  *currentSolution = arma::fft( *prevSolution );
  for ( unsigned int i=0;i<currentSolution->n_elem; i++ )
  {
    double kx = spatialFreq( i, currentSolution->n_elem );
    (*currentSolution)[i] *= kernel( kx );
  }
  *prevSolution = arma::ifft( *currentSolution );
}

void FFTSolver2D::refraction( unsigned int step )
{
  double z = guide->getZ( step ) + 0.5*stepZ;
  cdouble im(0.0,1.0);
  for ( unsigned int i=0;i<prevSolution->n_elem; i++ )
  {
    double delta, beta;
    double x = guide->getX(i);
    guide->getXrayMatProp( x, z, delta, beta );
    (*currentSolution)[i] = (*prevSolution)[i]*exp( -wavenumber*(beta+im*delta)*stepZ );
  }
}

double FFTSolver2D::spatialFreq( unsigned int indx, unsigned int size ) const
{
  if ( indx > size/2 )
  {
    indx = size-indx;
  }
  return 2.0*PI*indx/(stepX*size);
}

void FFTSolver2D::solveStep( unsigned int step )
{
  if ( callInit )
  {
    init();
    callInit = false;
  }
  propagate();
  refraction( step );
}

arma::cx_vec& FFTSolver2D::signalToFilter() const
{
  *filterSignal = currentSolution->subvec( 0, initialLength );
  return *filterSignal;
}
