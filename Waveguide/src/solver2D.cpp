#include "solver2D.hpp"
#include "waveGuideFDSimulation.hpp"
#include <complex>
#include <stdexcept>

using namespace std;

Solver2D::~Solver2D()
{
  if ( solution != NULL )
  {
    unsigned int Nx = guide->nodeNumberTransverse();
    unsigned int Nz = guide->nodeNumberLongitudinal();
    for ( unsigned int iz=0;iz<Nz;iz++ )
    {
      delete [] solution[iz];
    }
    delete [] solution;
  }
}

void Solver2D::setGuide( const WaveGuideFDSimulation &wg )
{
  guide = &wg;
  unsigned int Nx = guide->nodeNumberTransverse();
  unsigned int Nz = guide->nodeNumberLongitudinal();
  solution = new complex<double>*[Nz];
  for ( unsigned int iz=0;iz<Nz;iz++ )
  {
    solution[iz] = new complex<double>[Nx];
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
    solution[0][i] = values[i];
  }
}
