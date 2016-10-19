#include "solver1D.hpp"
#include <cassert>
#include <stdexcept>

using namespace std;

Solver1D::~Solver1D()
{
  delete solution;
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
  assert ( x >= x1 );
  assert ( x < x2 );
  unsigned int N = getEigenVectorSize();
  return N*( x- x1 )/( x2-x1 );
}

double Solver1D::getSolution( double x, unsigned int eigenmode ) const
{
  if (( x < x1) || ( x > x2 ))
  {
    return 0.0;
  }
  unsigned int indx = closestIndx(x);
  double dx = (x2-x1)/static_cast<double>(getEigenVectorSize());
  double x_low = x1+dx*indx;
  double x_high = x1+dx;
  return ((x-x_low)*getSolution()(indx,eigenmode) + (x_high-x)*getSolution()(indx+1,eigenmode))/dx;
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
  unsigned int ncol = solution->n_cols;
  solution->resize(ncol, nModes);
}
