#include "stdFDsolver.hpp"
#include <iostream>
#include "waveGuide.hpp"

using namespace std;

// Define the Lapack routine
extern "C" void dsteqr_( char *compz, int *N, double diag[], double subdiag[], double *Z,  int *LDZ, double *wrk, int *info);

void StandardFD::setLimits( double xlow, double xhigh )
{
  x1 = xlow;
  x2 = xhigh;
}

void StandardFD::fillJsonObj( Json::Value &obj ) const
{
  Solver1D::fillJsonObj( obj );
  obj["stepsize"] = stepsize;
  obj["xmin"] = x1;
  obj["xmax"] = x2;
}

void StandardFD::solve()
{
  int N = (x2-x1)/stepsize;
  double diag[N];
  double subdiag[N-1];

  // Fill matrix
  double hSq = stepsize*stepsize;
  for ( unsigned int i=0;i<N;i++ )
  {
    double x = x1 + stepsize*i;
    diag[i] = 2.0/hSq + waveguide->potential(x);
    if ( i < (N-1) )
    {
      subdiag[i] = -1.0/hSq;
    }
  }

  // Solve the symmetric tridiagonal system
  double eigvec[N*N];
  int info;
  char compz[] ="I";
  double wrk[2*N-2];
  dsteqr_(compz,&N,diag,subdiag,eigvec,&N,wrk, &info);
  if ( info == 0 )
  {
    clog << "System of equations solved successfully\n";
  }
  else if ( info < 0 )
  {
    cerr << "Error in argument " << -info << endl;
  }
  else
  {
    cerr << "Convergence problem. Return code " << info << endl;
  }

  // Copy the lowest eigenvalue
  eigenvalue = diag[0];
  for ( unsigned int i=0;i<N;i++ )
  {
    solution->push_back(eigvec[i]);
  }
}
