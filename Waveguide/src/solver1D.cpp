#include "solver1D.hpp"
#include <cassert>
#include <stdexcept>
#include <cmath>
#include "waveGuide.hpp"

using namespace std;

Solver1D::~Solver1D()
{
  delete solution;
  if ( solutionImag != NULL ) delete solutionImag;
}

void Solver1D::setLimits( double xlow, double xhigh )
{
  x1 = xlow;
  x2 = xhigh;
}

void Solver1D::fillJsonObj( Json::Value &obj ) const
{
  obj["name"] = name;
}

unsigned int Solver1D::closestIndx( double x ) const
{
  assert ( x > x1 );
  assert ( x < x2 );
  unsigned int N = getEigenVectorSize()-1;
  return N*( x - x1 )/( x2-x1 );
}

double Solver1D::getSolution( double x, unsigned int eigenmode ) const
{
  if (( x <= x1) || ( x >= x2 ))
  {
    return 0.0;
  }
  unsigned int indx = closestIndx(x);
  double dx = (x2-x1)/static_cast<double>(getEigenVectorSize());
  double x_low = x1+dx*indx;
  double x_high = x_low+dx;
  return ((x-x_low)*getSolution()(indx+1,eigenmode) + (x_high-x)*getSolution()(indx,eigenmode))/dx;
}

void Solver1D::addEigenmode( double eigenvalue, double eigvec[], unsigned int size )
{
  if ( nModes == 0 )
  {
    // This is the first eigenvector
    solution->set_size(size, maxNumberOfModes);
  }
  else if ( nModes >= maxNumberOfModes )
  {
    throw (runtime_error("Max number of modes reaches. Cannot load mode more"));
  }
  else if ( size != solution->n_rows )
  {
    throw( runtime_error("The dimensions of the eigenvector does not match the number of rows in the solution matrix."));
  }

  for ( unsigned int i=0;i<size;i++ )
  {
    (*solution)(i,nModes) = eigvec[i];
  }
  nModes++;
}

void Solver1D::loadingFinished()
{
  unsigned int nrows = solution->n_rows;
  solution->resize(nrows, nModes);
}

void Solver1D::printMode( unsigned int mode ) const
{
  unsigned int lineshiftEvery = 10;
  for ( unsigned int i=0;i<solution->n_rows;i++ )
  {
    cout << (*solution)(i,mode) << " ";
    if ( i%lineshiftEvery == 0 )
    {
      cout << "\n";
    }
  }
  cout << endl;
}

double Solver1D::normEigenvec( unsigned int mode ) const
{
  double norm = 0.0;
  // Trapezoidal integration
  norm += pow( (*solution)(0,mode), 2);
  norm += pow( (*solution)(solution->n_rows-1, mode), 2);
  for ( unsigned int i=1;i<solution->n_rows-1;i++ )
  {
    norm += 2.0*pow( (*solution)(i,mode), 2 );
  }
  return 0.5*norm*(x2-x1)/static_cast<double>(solution->n_rows);
}

void Solver1D::setGuide( const WaveGuide1DSimulation &guide )
{
  waveguide = &guide;
  if ( waveguide->useComplex() )
  {
    solutionImag = new arma::mat();
  }
}
