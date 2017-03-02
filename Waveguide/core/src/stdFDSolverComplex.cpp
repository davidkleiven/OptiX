#include "stdFDSolverComplex.hpp"
#include <complex>
#include <armadillo>
using namespace std;
typedef complex<double> cdouble;

extern "C" void zsteqr_( char *compz, int *N, cdouble *diag, cdouble *subdiag, cdouble *Z,  int *LDZ, cdouble *wrk, int *info);

void StdFDSolverComplex::solve()
{
  int N = (x2-x1)/stepsize;
  cx_mat A(N,N);
  cdouble diag[N];
  cdouble subdiag[N];
  // Fill matrix
  double hSq = stepsize*stepsize;
  for ( unsigned int i=0;i<N;i++ )
  {
    double x = x1 + stepsize*i;
    A(i,i) = 2.0/hSq + waveguide->potential(x);
    if ( i < (N-1) )
    {
      A(i,i+1) = -1.0/hSq;
      A(i+1,i) = A(i,i+1);
    }
  }

  // Solve the system
  cx_mat eigvec;
  cx_vec eigval;
  arma::eig_gen(eigval, eigvec, A);

  // Copy the solution
  solution->set_size(N,nModes);
  solutionImag->set_size(N,nModes);
  eigenvalues.clear();
  eigenvaluesImag.clear();
  for ( unsigned int i=0;i<nModes;i++ )
  {
    eigenvalues.push_back(eigval(i).real());
    eigenvaluesImag.push_back(eigval(i).imag());
    for ( unsigned int j=0;j<N;j++ )
    {
      (*solution)(j,i) = eigvec(j,i).real();
      (*solutionImag)(j,i) = eigvec(j,i).imag();
    }
  }
}
