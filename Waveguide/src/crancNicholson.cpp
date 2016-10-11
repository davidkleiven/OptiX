#include "crankNicholson.hpp"
#include "waveGuideFDSimulation.hpp"
#include <cassert>
#include <iostream>
#include <stdexcept>

using namespace std;

typedef complex<double> cdouble;

cdouble IMAG_UNIT(0.0,1.0);
CrankNicholson::~CrankNicholson(){};

void CrankNicholson::initValuesFromWaveGuide()
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
void CrankNicholson::solve()
{
  initValuesFromWaveGuide();
  for ( unsigned int iz=1;iz<Nz;iz++ )
  {
    solveCurrent( iz );
  }
}

void CrankNicholson::solveCurrent( unsigned int iz )
{
  assert( iz>=1 );
  cdouble *subdiag = new cdouble[Nx-1];
  cdouble *rhs = new cdouble[Nx];

  // Two useful dimensionless numbers
  double rho = stepZ/(wavenumber*stepX*stepX);
  double r = wavenumber*stepZ;

  cdouble alpha = 0.5*IMAG_UNIT/wavenumber;
  cdouble *diag = getSolution( iz ); // Pointer to where the next solution should be stored
  const cdouble *prevSol = getSolution( iz-1 );
  double z = zmin + iz*stepZ;
  for ( unsigned int ix=0;ix<Nx;ix++ )
  {
    double x = xmin + ix*stepX;
    double delta, beta;
    guide->getXrayMatProp( x, z, delta, beta );
    diag[ix] = 1.0 + 0.5*IMAG_UNIT*rho + beta*r + IMAG_UNIT*delta*r;

    if ( ix < Nx-1 )
    {
      subdiag[ix] = -0.25*IMAG_UNIT*rho;
    }

    // Fill right hand side
    if ( ix > 0 )
    {
      rhs[ix] = prevSol[ix-1];
    }

    if ( ix < Nx-1 )
    {
      rhs[ix] += prevSol[ix+1];
    }
    rhs[ix] *=  (0.25*IMAG_UNIT*rho);
    rhs[ix] -= 0.5*prevSol[ix]*IMAG_UNIT*rho;
    rhs[ix] += (1.0 - beta*r - IMAG_UNIT*delta*r)*prevSol[ix];
  }

  // Solve the tridiagonal system
  matrixSolver.solve( diag, subdiag, rhs, Nx);

  delete [] rhs;
  delete [] subdiag;
}

unsigned int CrankNicholson::rowColToIndx( unsigned int Nz, unsigned int ix, unsigned int iz )
{
  return ix*Nz+iz;
}

void CrankNicholson::indxToRowCol( unsigned int Nz, unsigned int indx, unsigned int &ix, unsigned int &iz )
{
  ix = indx/Nz;
  iz = indx%Nz;
}
