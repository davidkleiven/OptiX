#include "crankNicholson.hpp"
#include "waveGuideFDSimulation.hpp"
#include <cassert>
#include <iostream>
#include <stdexcept>
#include "paraxialEquation.hpp"
#include "borderTracker.hpp"
#include <cassert>
#include "boundaryCondition.hpp"

using namespace std;

typedef complex<double> cdouble;

cdouble IMAG_UNIT(0.0,1.0);
CrankNicholson::~CrankNicholson(){};

void CrankNicholson::solveStep( unsigned int iz )
{
  assert( iz>=1 );
  cdouble *subdiag = new cdouble[Nx-1];
  cdouble *rhs = new cdouble[Nx];
  cdouble *diag = new cdouble[Nx];

  // Two useful dimensionless numbers
  rho = stepZ/(wavenumber*stepX*stepX);
  double r = wavenumber*stepZ;

  double z = zmin + static_cast<double>(iz)*stepZ;
  BorderTracker* bTr = guide->getBorderTracker();

  if ( bTr != NULL )
  {
    bTr->locateBorder( z );
  }

  for ( unsigned int ix=0; ix<Nx; ix++ )
  {
    double x = guide->getX( ix );
    double xPrevShifted = x;
    double xShifted = x;
    if ( bTr != NULL )
    {
      xPrevShifted = bTr->getShiftedX(x, iz-1);
      xShifted = bTr->getShiftedX(x, iz);
    }

    double delta, beta;
    double deltaPrev, betaPrev;
    guide->getXrayMatProp( xShifted, z, delta, beta );
    guide->getXrayMatProp( xPrevShifted, z-stepZ, deltaPrev, betaPrev);

    Hpluss = eq->H(xShifted+0.5*stepX,z);
    Hminus = eq->H(xShifted-0.5*stepX,z);
    gval = eq->G(xShifted,z);
    double fval = eq->F(xShifted,z);
    double jval = eq->J(xShifted,z);
    double jvalPrev = eq->J(xPrevShifted,z-stepZ);
    delta -= jval;
    deltaPrev -= jvalPrev;

    diag[ix] = fval*1.0 + 0.25*( Hpluss+Hminus )*gval*IMAG_UNIT*rho + 0.5*(beta*r + IMAG_UNIT*delta*r);

    if ( ix < Nx-1 )
    {
      subdiag[ix] = -0.25*IMAG_UNIT*rho*Hminus*gval;
    }

    HplussPrev = eq->H(xPrevShifted+0.5*stepX,z-stepZ);
    HminusPrev = eq->H(xPrevShifted-0.5*stepX,z-stepZ);
    gvalPrev = eq->G(xPrevShifted,z-stepZ);
    double fvalPrev = eq->F(xPrevShifted,z-stepZ);

    // Fill right hand side
    cdouble left, center, right;
    if ( bTr != NULL )
    {
      bTr->threePointStencil( ix, iz, left, center, right);
    }
    else
    {
      if ( ix > 0 )
      {
        left = (*prevSolution)(ix-1);
      }
      else
      {
        left = 0.0;
      }
      center = (*prevSolution)(ix);
      if ( ix < Nx-1 )
      {
        right = (*prevSolution)(ix+1);
      }
      else
      {
        right = 0.0;
      }
    }

    rhs[ix] = left*HminusPrev*gvalPrev;
    rhs[ix] += right*HplussPrev*gvalPrev;

    rhs[ix] *=  (0.25*IMAG_UNIT*rho);
    rhs[ix] -= 0.25*center*IMAG_UNIT*rho*(HplussPrev+HminusPrev)*gval;
    rhs[ix] += (1.0*fvalPrev - 0.5*(betaPrev*r + IMAG_UNIT*deltaPrev*r) )*center;
  }

  applyBC( subdiag, diag, rhs );

  // Solve the tridiagonal system
  matrixSolver.solve( diag, subdiag, rhs, Nx);

  // Copy solution to matrix
  for ( unsigned int ix=0;ix<Nx;ix++ )
  {
    (*currentSolution)(ix) = diag[ix];
  }

  delete [] rhs;
  delete [] subdiag;
  delete [] diag;
}

void CrankNicholson::applyBC( cdouble subdiag[], cdouble diag[], cdouble rhs[] )
{
  switch ( boundaryCondition )
  {
    case BC_t::TRANSPARENT:
      if ( printBC ) clog << "Using transparent BC\n";
      applyTBC( subdiag, diag, rhs );
      break;
    default:
      if ( printBC ) clog << "Using default Dirichlet boundary condition\n";
  }
  printBC = false;
}

void CrankNicholson::applyTBC( cdouble subdiag[], cdouble diag[], cdouble rhs[] )
{
  const double ZERO = 1E-5;
  if ( abs( getLastSolution()(1) ) > ZERO )
  {
    applyTBCOneSide(subdiag, diag, rhs, 0, 1);
  }

  unsigned int N = getLastSolution().n_elem;
  if ( abs( getLastSolution()(N-2) ) > ZERO )
  {
    applyTBCOneSide( subdiag, diag, rhs, N-1, N-2 );
  }
}

void CrankNicholson::applyTBCOneSide( cdouble subdiag[], cdouble diag[], cdouble rhs[], unsigned int outer, unsigned int inner )
{
  cdouble im(0.0,1.0);
  cdouble ratio = getLastSolution()(outer)/getLastSolution()(inner);
  cdouble kdx = log( ratio )/im;
  if ( kdx.real() < 0.0 )
  {
    kdx.real(0.0);
  }

  double sign = guide->getX(outer) > guide->getX(inner) ? 1.0:-1.0;

  diag[outer] -= 0.25*im*rho*Hminus*gval*exp(im*kdx*sign);
  rhs[outer] += 0.25*im*rho*HminusPrev*gvalPrev*exp(im*kdx*sign)*getLastSolution()(outer);
}
