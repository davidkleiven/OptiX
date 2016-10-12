#include "solver2D.hpp"
#include "waveGuideFDSimulation.hpp"
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

void Solver2D::setGuide( const WaveGuideFDSimulation &wg )
{
  guide = &wg;
  unsigned int Nx = guide->nodeNumberTransverse();
  unsigned int Nz = guide->nodeNumberLongitudinal();
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
