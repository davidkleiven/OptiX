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
}

void Solver2D::setSimulator( ParaxialSimulation &wg )
{
  guide = &wg;
  unsigned int Nx = guide->nodeNumberTransverse();
  unsigned int Nz = guide->nodeNumberLongitudinal();
  if ( solution != NULL ) delete solution;
  solution = new arma::cx_mat(Nx,Nz);
}

void Solver2D::setLeftBC( const cdouble values[] )
{
  if ( solution == NULL )
  {
    throw (runtime_error("A solver must be given before setting the boundary conditions!"));
  }

  for ( unsigned int i=0;i<guide->nodeNumberTransverse();i++)
  {
    (*solution)(i,0) = values[i];
  }
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
