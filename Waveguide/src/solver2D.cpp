#include "solver2D.hpp"
#include "paraxialSimulation.hpp"
#include "paraxialEquation.hpp"
#include <complex>
#include <stdexcept>
#include <iostream>

using namespace std;

Solver2D::~Solver2D()
{
  if ( solution != NULL )
  {
    delete solution;
  }
  if ( prevSolution != NULL ) delete prevSolution;
  if ( currentSolution != NULL ) delete currentSolution;
}

void Solver2D::setSimulator( ParaxialSimulation &wg )
{
  Solver::setSimulator( wg );
  unsigned int Nx = guide->nodeNumberTransverse();
  unsigned int Nz = guide->nodeNumberLongitudinal();
  if ( solution != NULL ) delete solution;
  if ( prevSolution != NULL ) delete prevSolution;
  if ( currentSolution != NULL ) delete currentSolution;
  unsigned int downSampledNx = Nx/guide->transverseDiscretization().downsamplingRatio;
  solution = new arma::cx_mat(downSampledNx,Nz);
  prevSolution = new arma::cx_vec(Nx);
  currentSolution = new arma::cx_vec(Nx);

  if ( downSampledNx != Nx )
  {
    filter.setSourceSize( Nx );
    filter.setTargetSize( downSampledNx );
    filter.computeFilterCoefficients( kernel );
  }
}

void Solver2D::setLeftBC( const cdouble values[] )
{
  if ( solution == NULL )
  {
    throw (runtime_error("A solver must be given before setting the boundary conditions!"));
  }

  for ( unsigned int i=0;i<guide->nodeNumberTransverse();i++)
  {
    (*currentSolution)(i) = values[i];
  }
  copyCurrentSolution(0);
}

void Solver2D::setXBC( const cdouble valuesTop[], const cdouble valuesBottom[] )
{
  if ( solution == NULL )
  {
    throw (runtime_error("A solver must be given before setting the boundary conditions!"));
  }

  unsigned int row = solution->n_rows-1;
  for ( unsigned int i=0;i<guide->nodeNumberLongitudinal();i++ )
  {
    (*solution)(row,i) = valuesTop[i];
    (*solution)(0,i) = valuesBottom[i];
  }
}

void Solver2D::realPart( double *realsol ) const
{
  realOrImagPart(realsol, Comp_t::REAL);
}
void Solver2D::imagPart( double *imagsol ) const
{
  realOrImagPart(imagsol, Comp_t::IMAG);
}

void Solver2D::realOrImagPart( double *compsolution, Comp_t comp ) const
{
  if ( solution == NULL )
  {
    throw (runtime_error("No equation system has been solved!"));
  }

  unsigned int Nx = guide->nodeNumberTransverse();
  unsigned int Nz = guide->nodeNumberLongitudinal();
  for ( unsigned int ix=0; ix<Nx; ix++ )
  {
    for ( unsigned int iz=0; iz < Nz; iz++ )
    {
      switch ( comp )
      {
        case Comp_t::REAL:
          compsolution[ix*Nz+iz] = (*solution)(ix,iz).real();
          break;
        case Comp_t::IMAG:
          compsolution[ix*Nz+iz] = (*solution)(ix,iz).imag();
          break;
      }
    }
  }
}

void Solver2D::fillInfo( Json::Value &obj ) const
{
  obj["name"] = name;
}

bool Solver2D::importHDF5( const string &fname )
{
  return solution->load( fname.c_str(), arma::hdf5_binary );
}

bool Solver2D::importHDF5( const string &realpart, const string &amplitude )
{
  // Load the real part
  bool status = solution->load( realpart.c_str(), arma::hdf5_binary );
  if ( !status ) return status;

  arma::mat amp;
  status = amp.load( amplitude.c_str(), arma::hdf5_binary);
  if ( !status ) return status;

  clog << "Warning: Currently there is a bug such that this function only populates the real part of the solution matrix!\n";
  // The imaginary part can now have a sign error. Should not matter
  // TODO: The next line fails. Why?
  //solution->set_imag( arma::sqrt(arma::pow(amp,2) - arma::pow(arma::real(*solution),2)) );
  return status;
}

void Solver2D::getField(arma::mat &field ) const
{
  field.set_size( solution->n_rows, solution->n_cols );
  double k = guide->getWavenumber();
  for ( unsigned int iz=0;iz<guide->nodeNumberLongitudinal();iz++ )
  {
    for ( unsigned int ix=0;ix<guide->nodeNumberTransverse();ix++ )
    {
      field(ix,iz) = (*solution)(ix,iz).real();//*eq->phaseFactor(k, guide->getZ(iz)) ).real();
    }
  }
}

void Solver2D::getPhase( arma::mat &phase ) const
{
  phase.set_size( solution->n_rows, solution->n_cols );
  for ( unsigned int iz=0;iz<guide->nodeNumberLongitudinal(); iz++ )
  {
    for ( unsigned int ix=0;ix<guide->nodeNumberTransverse(); ix++ )
    {
      phase(ix,iz) = arg( (*solution)(ix,iz) );
    }
  }
}

void Solver2D::copyCurrentSolution( unsigned int step )
{
  // Copy results to the previous array
  for ( unsigned int i=0;i<currentSolution->n_elem;i++ )
  {
    (*prevSolution)(i) = (*currentSolution)(i);
  }

  if ( guide->transverseDiscretization().downsamplingRatio > 1 )
  {
      visa::ArmaGetter<cdouble,visa::ArmaMatrix_t::VECTOR> getter;
      filter.filterArray( signalToFilter(), getter );
  }

  // Downsample array
  double delta = static_cast<double>( signalToFilter().n_elem )/static_cast<double>( solution->n_rows );
  for ( unsigned int i=0;i<solution->n_rows;i++ )
  {
    (*solution)(i,step) = signalToFilter()(delta*i);
  }
}

void Solver2D::filterInLongitudinalDirection()
{
  if ( guide->longitudinalDiscretization().downsamplingRatio == 1 ) return;

  unsigned int downSampledNz = solution->n_cols/guide->longitudinalDiscretization().downsamplingRatio;
  filter.setTargetSize( downSampledNz );
  filter.setSourceSize( solution->n_cols );
  filter.computeFilterCoefficients( kernel );

  visa::ArmaGetter<cdouble, visa::ArmaMatrix_t::ROW> getter;
  // Filter the matrix
  for ( unsigned int i=0;i<solution->n_rows;i++ )
  {
    filter.filterArray( *solution, getter );
  }

  arma::cx_mat *downSampledMatrix = new arma::cx_mat( solution->n_rows, downSampledNz );
  double delta = static_cast<double>( solution->n_cols )/static_cast<double>( downSampledNz );
  //double delta = static_cast<double>( solution->n_rows )/static_cast<double>( downSampledNz );

  for ( unsigned int i=0;i<downSampledNz;i++ )
  {
    for ( unsigned int j=0;j<solution->n_rows;j++ )
    {
      (*downSampledMatrix)(j,i) = (*solution)(j, delta*i);
    }
  }
}

void Solver2D::solve()
{
  if ( eq == NULL )
  {
    throw ( runtime_error("No paraxial equation object given!") );
  }

  if ( bc == NULL )
  {
    throw ( runtime_error("No boundary conditions given!") );
  }

  // Assert that the solution matrix is allocated
  assert( solution != NULL );
  assert( guide != NULL );

  initValuesFromWaveGuide();
  for ( unsigned int iz=1;iz<Nz;iz++ )
  {
    step();
  }
  filterInLongitudinalDirection();
  downSampleLongitudinalDirection();
}

void Solver2D::step()
{
  if ( currentStep == 1 ) initValuesFromWaveGuide();

  solveStep( currentStep );
  copyCurrentSolution( currentStep );
  currentStep++;
}

void Solver2D::initValuesFromWaveGuide()
{
  if ( guide == NULL )
  {
    throw( runtime_error("No waveguide specified!"));
  }
  Nx = guide->nodeNumberTransverse();
  Nz = guide->nodeNumberLongitudinal();
  stepX = guide->transverseDiscretization().step;
  xmin = guide->transverseDiscretization().min;
  stepZ = guide->longitudinalDiscretization().step;
  zmin = guide->longitudinalDiscretization().min;
  wavenumber = guide->getWavenumber();
}

void Solver2D::downSampleLongitudinalDirection()
{
  unsigned int Nz = guide->nodeNumberLongitudinal()/guide->longitudinalDiscretization().downsamplingRatio;
  arma::cx_mat *copy = new arma::cx_mat( solution->n_rows, Nz );
  double delta = static_cast<double>( solution->n_cols )/static_cast<double>( Nz );
  for ( unsigned int iz=0;iz<Nz;iz++ )
  {
    unsigned int indx = iz*delta+delta/2.0;
    for ( unsigned int ix=0;ix<copy->n_rows;ix++ )
    {
      (*copy)(ix,iz) = (*solution)(ix,indx);
    }
  }
  delete solution;
  solution = copy;
}
