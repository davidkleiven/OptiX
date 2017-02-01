#include "postProcessMod.hpp"
#include "paraxialSimulation.hpp"
#include "solver.hpp"
#include <cmath>
#include <fftw3.h>
#include <armadillo>
#include <gsl/gsl_fit.h>
#include <complex>

using namespace std;
typedef complex<double> cdouble;

const double PI = acos(-1.0);
void post::FarField::setAngleRange( double angMin, double angMax )
{
  phiMin = angMin;
  phiMax = angMax;
}

void post::FarField::result( const Solver &solver, arma::vec &res )
{
  // Extract the last column of the solution matrix
  arma::cx_vec exitField = solver.getLastSolution();
  arma::cx_vec paddedSignal;
  double zpos = sim->longitudinalDiscretization().max;

  if ( signalLength < exitField.n_elem )
  {
    paddedSignal = exitField;
  }
  else
  {
    paddedSignal.set_size(signalLength);
    paddedSignal.fill(0.0);

    unsigned int start = signalLength/2 - exitField.n_elem/2;
    // Pad the signal
    for ( int i=0;i<paddedSignal.n_elem;i++ )
    {
      int indx = i-static_cast<int>(start);
      double x = sim->getX( indx );
      paddedSignal[i] = sim->padExitField( x, zpos );
    }
    // Fill in the signal on the center
    for ( unsigned int i=0;i<exitField.n_elem;i++ )
    {
      paddedSignal(start+i) = exitField(i);
    }
  }

  arma::cx_vec ft = arma::fft( paddedSignal );

  res = arma::abs(ft)/sqrt(ft.n_elem);

  // Shift the FFT
  fftshift( res );

  // Reduce the result array
  reduceArray( res );
}

template<class T>
void post::FarField::fftshift( arma::Col<T> &array )
{
  unsigned int N = array.n_elem;
  for ( unsigned int i=0;i<N/2;i++ )
  {
    T copy = array(i);
    array(i) = array(i+N/2);
    array(i+N/2) = copy;
  }
}
template void post::FarField::fftshift<double>( arma::Col<double> &array );
template void post::FarField::fftshift<cdouble>( arma::Col<cdouble> &array );

unsigned int post::FarField::farFieldAngleToIndx( double angle, unsigned int size ) const
{
  double angleRad = angle*PI/180.0;
  int n = size*sim->getWavenumber()*sim->transverseDiscretization().step*sin(angleRad)/(2.0*PI);

  //int n = size*sim->getWavenumber()*sin(angleRad)/(2.0*PI*sim->transverseDiscretization().step);
  if ( abs(n) >= size/2 )
  {
    if ( angle > 0.0 ) return size-1;
    return 0;
  }
  return n+static_cast<int>(size)/2;
}

template<class T>
void post::FarField::reduceArray( arma::Col<T> &res ) const
{
  unsigned int indxMin = farFieldAngleToIndx( phiMin, res.n_elem );
  unsigned int indxMax = farFieldAngleToIndx( phiMax, res.n_elem );
  res = res.subvec( indxMin, indxMax );
}

template void post::FarField::reduceArray<double>( arma::Col<double> &res ) const;
template void post::FarField::reduceArray<cdouble>( arma::Col<cdouble> &res ) const;

void post::FarField::reduceArray( arma::mat &res ) const
{
  unsigned int indxMin = farFieldAngleToIndx( phiMin, res.n_rows );
  unsigned int indxMax = farFieldAngleToIndx( phiMax, res.n_rows );
  clog << indxMin << " " << indxMax << endl;
  res = res.submat( indxMin, indxMin, indxMax, indxMax );
}

void post::FarField::addAttrib( vector<H5Attr> &attr ) const
{
  attr.push_back( makeAttr("phiMin", phiMin) );
  attr.push_back( makeAttr("phiMax", phiMax) );
  double qmax = phiMax*sim->getWavenumber()*PI/180.0;
  double qmin = phiMin*sim->getWavenumber()*PI/180.0;
  attr.push_back( makeAttr("qmin", qmin) );
  attr.push_back( makeAttr("qmax", qmax) );
}

void post::FarField::result( const Solver &solver, arma::mat &res )
{
  unsigned int Nx = signalLength < res.n_rows ? res.n_rows:signalLength;
  unsigned int Ny = signalLength < res.n_cols ? res.n_cols:signalLength;

  // Perform FFT over columns
  arma::cx_vec pad( Nx );
  arma::cx_vec ft( Nx );
  pad.fill(0.0);
  unsigned int indxMin = farFieldAngleToIndx( phiMin, pad.n_elem );
  unsigned int indxMax = farFieldAngleToIndx( phiMax, pad.n_elem );
  assert( indxMax >= indxMin );

  unsigned int size = indxMax-indxMin+1;
  clog << "Size: " << size << endl;
  if ( size == 0 )
  {
    cout << "The requested far field size is zero!\n";
    return;
  }
  res.set_size( size, size );

  arma::cx_mat solution = solver.getLastSolution3D();
  if ( reference != NULL )
  {
    solution -= *reference;
  }

  arma::cx_mat temporary(size, solution.n_cols );
  temporary.fill( 0.0 );

  fftw_complex* data = reinterpret_cast<fftw_complex*>( pad.memptr() );

  // TODO: If mysterious seg. faults occure, check if the pad signal needs to have to extra elements when performing in-place transform
  fftw_plan plan = fftw_plan_dft_1d( pad.n_elem, data, data, FFTW_FORWARD, FFTW_ESTIMATE );

  // FFT over columns
  for ( unsigned int i=0;i<solution.n_cols;i++ )
  {
    pad.subvec( 0, solution.n_rows-1 ) = solution.col(i);
    padSignal( pad );
    fftw_execute( plan );
    ft = pad; // Perform copy to avoid changing the FFTW plan
    fftshift( ft );
    reduceArray( ft );
    //temporary.insert_cols(i, ft);
    temporary.col(i) = ft;
    pad.fill(0.0);
  }

  // FFT over rows
  assert( Nx == Ny );
  pad.set_size( Ny );
  pad.fill( 0.0 );
  for ( unsigned int i=0;i<temporary.n_rows;i++ )
  {
    pad.subvec( 0, temporary.n_cols-1 ) = temporary.row(i).t();
    padSignal( pad );
    fftw_execute( plan );
    ft = pad; // Copy to not change the FFTW plan
    fftshift( ft );
    reduceArray( ft );
    //res.insert_rows( i, arma::pow( arma::abs( ft ), 2 ).t() );
    res.row(i) = arma::pow( arma::abs( ft ), 2 ).t()/signalLength;
    pad.fill(0.0);
  }
  fftw_destroy_plan( plan );
}

void post::FarField::padSignal( arma::cx_vec &zeroPadded ) const
{
  switch ( padding )
  {
    case Pad_t::ZERO:
      return;
    case Pad_t::LINEAR_FIT:
      clog << "Warning! Linear padding has not been tested properly and may not work!\n";
      linearPadding( zeroPadded );
      return;
    case Pad_t::CONSTANT:
      constantPadding( zeroPadded );
      return;
  }
}

void post::FarField::linearPadding( arma::cx_vec &zeroPadded ) const
{
  unsigned int nInFit = 10;
  // TODO: Is this safe?
  double *dataStart = reinterpret_cast<double*>( zeroPadded.memptr() );
  size_t stride = 2; // Skip imaginary parts
  double xmin = sim->getX(0);
  double xmax = sim->getX(nInFit);
  arma::vec x = arma::linspace(0.0,xmax-xmin,nInFit);
  double c0, c1, cov00, cov01, cov11, sumsq;
  gsl_fit_linear( x.memptr(), 1, dataStart, stride, nInFit, &c0, &c1, &cov00, &cov01, &cov11, &sumsq );

  dataStart += 1;
  double d0, d1;
  gsl_fit_linear( x.memptr(), 1, dataStart, stride, nInFit, &d0, &d1, &cov00, &cov01, &cov11, &sumsq );

  // Pad from beginning
  unsigned int N = sim->nodeNumberTransverse();
  assert( N <= zeroPadded.n_elem );
  unsigned int added = zeroPadded.n_elem - N;
  unsigned int padStart = zeroPadded.n_elem - added/2;
  cdouble a0 = c0+1i*d0;
  cdouble a1 = c1+1i*d1;
  double dx = x[1] - x[0];
  for ( int i=padStart;i<zeroPadded.n_elem;i++ )
  {
    double currX = x[0] + dx*(i-zeroPadded.n_elem+1);
    zeroPadded[i] = a0 + a1*currX;
  }

  // Pad from end
  xmin = sim->getX( sim->nodeNumberTransverse()-nInFit );
  xmax = sim->getX( sim->nodeNumberTransverse() );
  x = arma::linspace( xmin, xmax, nInFit );
  dataStart = reinterpret_cast<double*>( zeroPadded.memptr() );
  gsl_fit_linear( x.memptr(), 1, dataStart, stride, nInFit, &c0, &c1, &cov00, &cov01, &cov11, &sumsq );

  dataStart += 1;
  gsl_fit_linear( x.memptr(), 1, dataStart, stride, nInFit, &d0, &d1, &cov00, &cov01, &cov11, &sumsq );
  padStart = sim->nodeNumberTransverse();
  a0 = c0+1i*d0;
  a1 = c1+1i*d1;
  for ( int i=padStart;i<padStart+added/2;i++)
  {
    double currX = x[x.n_elem-1] + (i-padStart)*dx;
    zeroPadded[i] = a0 + a1*currX;
  }
}

void post::FarField::constantPadding( arma::cx_vec &zeroPadded ) const
{
  unsigned int N = sim->nodeNumberTransverse();
  for ( unsigned int i=N;i<zeroPadded.n_elem;i++ )
  {
    zeroPadded[i] = zeroPadded[0];
  }
}
