#include "crankNicholson.hpp"
#include "waveGuideFDSimulation.hpp"
#include <cassert>
#include <iostream>
#include <stdexcept>
#include "paraxialEquation.hpp"

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
  if ( eq == NULL )
  {
    throw ( runtime_error("No paraxial equation object given!") );
  }

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
    double deltaPrev, betaPrev;
    guide->getXrayMatProp( x, z, delta, beta );
    guide->getXrayMatProp( x, z-stepZ, deltaPrev, betaPrev);

    double Hpluss = eq->H(x+0.5*x,z);
    double Hminus = eq->H(x-0.5*x,z);
    double gval = eq->G(x,z);
    double fval = eq->F(x,z);
    double jval = eq->J(x,z);
    double jvalPrev = eq->J(x,z-stepZ);
    delta -= jval;
    deltaPrev -= jvalPrev;

    diag[ix] = fval*1.0 + 0.25*( Hpluss+Hminus )*gval*IMAG_UNIT*rho + 0.5*(beta*r + IMAG_UNIT*delta*r);

    if ( ix < Nx-1 )
    {
      subdiag[ix] = -0.25*IMAG_UNIT*rho*Hminus*gval;
    }

    double HplussPrev = eq->H(x+0.5*x,z-stepZ);
    double HminusPrev = eq->H(x-0.5*x,z-stepZ);
    double gvalPrev = eq->G(x,z-stepZ);
    double fvalPrev = eq->F(x,z-stepZ);

    // Fill right hand side
    if ( ix > 0 )
    {
      rhs[ix] = (*solution)(ix-1,iz-1)*Hminus*gval;
    }
    else
    {
      rhs[ix] = guide->transverseBC(z, WaveGuideFDSimulation::Boundary_t::BOTTOM)*(Hminus*gval+HminusPrev*gvalPrev); // Make sure that it is not a random value
    }

    if ( ix < Nx-1 )
    {
      rhs[ix] += (*solution)(ix+1,iz-1)*HplussPrev*gvalPrev;
    }
    else
    {
      rhs[ix] += guide->transverseBC(z, WaveGuideFDSimulation::Boundary_t::TOP)*(Hpluss*gval + HplussPrev*gvalPrev);
    }

    rhs[ix] *=  (0.25*IMAG_UNIT*rho);
    rhs[ix] -= 0.25*(*solution)(ix,iz-1)*IMAG_UNIT*rho*(HplussPrev+HminusPrev)*gval;
    rhs[ix] += (1.0*fvalPrev - 0.5*(betaPrev*r + IMAG_UNIT*deltaPrev*r) )*(*solution)(ix,iz-1);
  }

  // Solve the tridiagonal system
  matrixSolver.solve( diag, subdiag, rhs, Nx);

  // Copy solution to matrix
  for ( unsigned int ix=0;ix<Nx;ix++ )
  {
    (*solution)(ix,iz) = diag[ix];
  }

  delete [] rhs;
  delete [] subdiag;
  delete [] diag;
}
