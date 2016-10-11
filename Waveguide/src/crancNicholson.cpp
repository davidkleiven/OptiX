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
  cdouble *diag = new cdouble[Nx];

  // Two useful dimensionless numbers
  double rho = stepZ/(wavenumber*stepX*stepX);
  double r = wavenumber*stepZ;

  double z = zmin + static_cast<double>(iz)*stepZ;
  for ( unsigned int ix=0; ix<Nx; ix++ )
  {
    double x = xmin + static_cast<double>(ix)*stepX;
    double delta, beta;
    guide->getXrayMatProp( x, z, delta, beta );
    diag[ix] = 1.0 + 0.5*IMAG_UNIT*rho + 0.5*(beta*r + IMAG_UNIT*delta*r);

    if ( ix < Nx-1 )
    {
      subdiag[ix] = -0.25*IMAG_UNIT*rho;
    }

    // Fill right hand side
    if ( ix > 0 )
    {
      rhs[ix] = solution[ix-1][iz-1];
    }
    else
    {
      rhs[ix] = 0.0; // Make sure that it is not a random value
    }

    if ( ix < Nx-1 )
    {
      rhs[ix] += solution[ix+1][iz-1];
    }
    rhs[ix] *=  (0.25*IMAG_UNIT*rho);
    rhs[ix] -= 0.5*solution[ix][iz-1]*IMAG_UNIT*rho;
    rhs[ix] += (1.0 - 0.5*(beta*r + IMAG_UNIT*delta*r) )*solution[ix][iz-1];
  }

  // Solve the tridiagonal system
  matrixSolver.solve( diag, subdiag, rhs, Nx);

  // Copy solution to matrix
  for ( unsigned int ix=0;ix<Nx;ix++ )
  {
    solution[ix][iz] = diag[ix];
  }

  delete [] rhs;
  delete [] subdiag;
  delete [] diag;
}
