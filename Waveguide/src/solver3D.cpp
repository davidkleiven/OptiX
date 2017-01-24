#include "solver3D.hpp"
#include <cassert>
#include "paraxialSimulation.hpp"

using namespace std;

void Solver3D::setSimulator( ParaxialSimulation &sim )
{
  Solver::setSimulator( sim );
  Nx = guide->nodeNumberTransverse();
  Nz = guide->nodeNumberLongitudinal();
  Ny = Nx;

  unsigned int downSampledX = Nx/guide->transverseDiscretization().downsamplingRatio;
  unsigned int downSampledZ = Nz/guide->longitudinalDiscretization().downsamplingRatio;

  // Deallocate if already allocated
  if ( solution != NULL ) delete solution;
  if ( prevSolution != NULL ) delete prevSolution;
  if ( currentSolution != NULL ) delete currentSolution;

  prevSolution = new arma::cx_mat(Nx,Ny);
  currentSolution = new arma::cx_mat(Nx, Ny);
  solution = new arma::cx_cube( Nx, Ny, Nz );

  if ( downSampledX != Nx )
  {
    filter.setSourceSize( Nx );
    filter.setTargetSize( downSampledX );
    filter.computeFilterCoefficients( kernel );
  }
}

void Solver3D::filterTransverse( arma::cx_mat &mat )
{
  assert( Nx == Ny ); // If this is the case, the filter coefficients does not need to be recomputed

  visa::ArmaGetter<cdouble, visa::ArmaMatrix_t::COL> colGetter;
  for ( unsigned int i=0;i<mat.n_rows;i++ )
  {
    colGetter.fixedIndx = i;
    filter.filterArray( mat, colGetter );
  }

  visa::ArmaGetter<cdouble, visa::ArmaMatrix_t::ROW> rowGetter;
  for ( unsigned int i=0;i<mat.n_cols;i++ )
  {
    rowGetter.fixedIndx = i;
    filter.filterArray( mat, rowGetter );
  }
}

void Solver3D::copyCurrentSolution( unsigned int step )
{
  assert( currentSolution->n_rows == prevSolution->n_rows );
  assert( currentSolution->n_cols == prevSolution->n_cols );

  for ( unsigned int i=0;i<currentSolution->n_rows; i++ )
  {
    for ( unsigned int j=0;j<currentSolution->n_cols;j++ )
    {
      (*prevSolution)(i,j) = (*currentSolution)(i,j);
    }
  }

  if (( currentSolution->n_rows != solution->n_rows ) || ( currentSolution->n_cols != solution->n_cols ))
  {
    filterTransverse( *currentSolution );
  }

  double deltaX = static_cast<double>( currentSolution->n_rows )/static_cast<double>( solution->n_rows );
  double deltaY = static_cast<double>( currentSolution->n_cols )/static_cast<double>( solution->n_cols );

  unsigned int maxX = deltaX*( solution->n_rows - 1 );
  unsigned int maxY = deltaY*( solution->n_cols - 1 );
  assert( maxX < currentSolution->n_rows );
  assert( maxY < currentSolution->n_cols );

  // Downsample the array
  for ( unsigned int i=0;i<solution->n_rows;i++ )
  {
    for ( unsigned int j=0;j<solution->n_cols;j++ )
    {
      (*solution)(i,j,step) = (*currentSolution)( i*deltaX, j*deltaY );
    }
  }
}
