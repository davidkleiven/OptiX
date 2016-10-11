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
    unsigned int Nx = guide->nodeNumberTransverse();
    for ( unsigned int ix=0;ix<Nx;ix++ )
    {
      delete [] solution[ix];
    }
    delete [] solution;
  }
}

void Solver2D::setGuide( const WaveGuideFDSimulation &wg )
{
  guide = &wg;
  unsigned int Nx = guide->nodeNumberTransverse();
  unsigned int Nz = guide->nodeNumberLongitudinal();
  solution = new complex<double>*[Nx];
  for ( unsigned int ix=0;ix<Nx;ix++ )
  {
    solution[ix] = new complex<double>[Nz];
  }
}

void Solver2D::setLeftBC( const cdouble values[] )
{
  if ( solution == NULL )
  {
    throw (runtime_error("A solver must be given before setting the boundary conditions!"));
  }
  for ( unsigned int i=0;i<guide->nodeNumberTransverse();i++ )
  {
    solution[i][0] = values[i];
  }
}

void Solver2D::realPart( double **realsol ) const
{
  realOrImagPart(realsol, Comp_t::REAL);
}
void Solver2D::imagPart( double **imagsol ) const
{
  realOrImagPart(imagsol, Comp_t::IMAG);
}

void Solver2D::realOrImagPart( double **compsolution, Comp_t comp ) const
{
  if ( solution == NULL )
  {
    throw (runtime_error("No equation system has been solved!"));
  }

  for ( unsigned int ix=0; ix<guide->nodeNumberTransverse(); ix++ )
  {
    for ( unsigned int iz=0; iz < guide->nodeNumberLongitudinal(); iz++ )
    {
      switch ( comp )
      {
        case Comp_t::REAL:
          compsolution[ix][iz] = solution[ix][iz].real();
          break;
        case Comp_t::IMAG:
          compsolution[ix][iz] = solution[ix][iz].imag();
          break;
      }
    }
  }
}

void Solver2D::fillInfo( Json::Value &obj ) const
{
  obj["name"] = name;
}
