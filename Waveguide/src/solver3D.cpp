#include "solver3D.hpp"
#include <cassert>
#include "paraxialSimulation.hpp"
#include <omp.h>
#include <ctime>
#include <chrono>
#include <iostream>

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

  #pragma omp parallel for
  for ( unsigned int i=0;i<mat.n_rows;i++ )
  {
    colGetter.fixedIndx = i;
    filter.filterArray( mat, colGetter );
  }

  unsigned int maxThreads = omp_get_num_threads();
  vector< visa::ArmaGetter<cdouble, visa::ArmaMatrix_t::ROW> > rowGetters;
  for ( unsigned int i=0;i<maxThreads;i++ )
  {
    rowGetters.push_back( visa::ArmaGetter<cdouble, visa::ArmaMatrix_t::ROW>() );
  }

  #pragma omp parallel
  {
    unsigned int id = omp_get_thread_num();
    #pragma omp for
    for ( unsigned int i=0;i<mat.n_cols;i++ )
    {
      rowGetters[id].fixedIndx = i;
      filter.filterArray( mat, rowGetters[id] );
    }
  }
}

void Solver3D::copyCurrentSolution( unsigned int step )
{
  assert( currentSolution->n_rows == prevSolution->n_rows );
  assert( currentSolution->n_cols == prevSolution->n_cols );

  *prevSolution = *currentSolution;

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
  #pragma omp parallel for
  for ( unsigned int i=0;i<solution->n_rows*solution->n_cols;i++ )
  {
    unsigned int row = i%solution->n_rows;
    unsigned int col = i/solution->n_rows;
    (*solution)(row,col,step) = (*currentSolution)( row*deltaX, col*deltaY );
  }

}

void Solver3D::setInitialConditions( const arma::cx_mat &values )
{

  if (( currentSolution == NULL ) || ( prevSolution == NULL ) || ( solution == NULL ))
  {
    throw ( runtime_error("The function setSimulator needs to be called before setInitialConditions!") );
  }

  if (( values.n_rows != currentSolution->n_rows ) || ( values.n_cols != currentSolution->n_cols ))
  {
    stringstream ss;
    ss << "Dimension of matrices does not match!\n";
    ss << "Given: Nrows: " << values.n_rows << " Ncols: " << values.n_cols;
    ss << "\nRequired: Nrows: " << currentSolution->n_rows << " Ncols: " << currentSolution->n_cols;
    throw( runtime_error( ss.str() ) );
  }

  (*currentSolution) = values;
  copyCurrentSolution( 0 );
}

void Solver3D::step()
{
  solveStep( currentStep );
  copyCurrentSolution( currentStep++ );
}

void Solver3D::solve()
{
  auto lastTime = chrono::steady_clock::now();
  auto now = lastTime;
  for ( unsigned int i=1;i<guide->nodeNumberLongitudinal();i++ )
  {
    step();

    now = chrono::steady_clock::now();
    chrono::duration<double> elapsedSec = now - lastTime;
    if ( elapsedSec > chrono::duration<double>(secBetweenStatusMessage) )
    {
      clog << "Propagation step: " << i << " of " << guide->nodeNumberLongitudinal() << endl;
      lastTime = chrono::steady_clock::now();
    }
  }

  // TODO: Should one support downsampling in 3D. This will require an extra copy
  //filterInLongitudinalDirection();
  //downSampleLongitudinalDirection();
}

const arma::cx_cube& Solver3D::getSolution3D() const
{
  return *solution;
}
