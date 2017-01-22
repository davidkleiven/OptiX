#include "fftSolver2D.hpp"
#include <cmath>
#include "paraxialSimulation.hpp"
#include <iostream>

using namespace std;
const double PI = acos(-1.0);

FFTSolver2D::~FFTSolver2D()
{
  fftw_destroy_plan( ftforw );
  fftw_destroy_plan( ftback );
}

cdouble FFTSolver2D::kernel( double kx ) const
{
  cdouble A( 0.0, 0.5*stepZ/wavenumber );
  return exp( -A*kx*kx );
}

void FFTSolver2D::propagate()
{
  //*currentSolution = arma::fft( *prevSolution );
  //clog << arma::max( arma::abs( *prevSolution) ) << endl;
  fftw_execute( ftforw ); // FFT(prevSolution) -> currentSolution
  for ( unsigned int i=0;i<prevSolution->n_elem; i++ )
  {
    double kx = spatialFreq( i, prevSolution->n_elem );
    //(*prevSolution)[i] *= kernel( kx );
    (*currentSolution)[i] *= kernel( kx );
  }
  //*prevSolution = arma::ifft( *currentSolution );
  fftw_execute( ftback ); // IFFT(currentSolution) -> prevSolution
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

    // FFTW3: Divide by length to normalize
    double normalization = prevSolution->n_elem;
    //normalization = 1.0;
    (*currentSolution)[i] = (*prevSolution)[i]*exp( -wavenumber*(beta+im*delta)*stepZ )/normalization;
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
  //planInitialized = true;
  if ( !planInitialized )
  {
    prev = reinterpret_cast<fftw_complex*>( prevSolution->memptr() );
    curr = reinterpret_cast<fftw_complex*>( currentSolution->memptr() );
    ftforw = fftw_plan_dft_1d( prevSolution->n_elem, prev, curr, FFTW_FORWARD, FFTW_ESTIMATE );
    ftback = fftw_plan_dft_1d( prevSolution->n_elem, curr, prev, FFTW_BACKWARD, FFTW_ESTIMATE );
    planInitialized = true;
  }
  propagate();
  refraction( step );
}
