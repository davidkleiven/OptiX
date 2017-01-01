#include "postProcessMod.hpp"
#include "paraxialSimulation.hpp"
#include "solver2D.hpp"
#include <cmath>

using namespace std;

const double PI = acos(-1.0);
void post::FarField::setAngleRange( double angMin, double angMax )
{
  phiMin = angMin;
  phiMax = angMax;
}

void post::FarField::result( const Solver2D &solver, arma::vec &res )
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
  unsigned int N = res.n_elem;
  for ( unsigned int i=0;i<N/2;i++ )
  {
    double copy = res(i);
    res(i) = res(i+N/2);
    res(i+N/2) = copy;
  }

  // Reduce the result array
  reduceArray( res );
}

unsigned int post::FarField::farFieldAngleToIndx( double angle, const arma::vec &res ) const
{
  double angleRad = angle*PI/180.0;
  int n = res.n_elem*sim->getWavenumber()*sim->transverseDiscretization().step*sin(angleRad)/(2.0*PI);
  if ( abs(n) >= res.n_elem/2 )
  {
    if ( angle > 0.0 ) return res.n_elem-1;
    return 0;
  }
  return n+static_cast<int>(res.n_elem)/2;
}

void post::FarField::reduceArray( arma::vec &res ) const
{
  unsigned int indxMin = farFieldAngleToIndx( phiMin, res );
  unsigned int indxMax = farFieldAngleToIndx( phiMax, res );
  res = res.subvec( indxMin, indxMax );
}

void post::FarField::addAttrib( vector<H5Attr> &attr ) const
{
  attr.push_back( makeAttr("phiMin", phiMin) );
  attr.push_back( makeAttr("phiMax", phiMax) );
}
