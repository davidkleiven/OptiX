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
  double rho = stepZ/(wavenumber*stepX*stepX);
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

    double Hpluss = eq->H(xShifted+0.5*stepX,z);
    double Hminus = eq->H(xShifted-0.5*stepX,z);
    double gval = eq->G(xShifted,z);
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

    double HplussPrev = eq->H(xPrevShifted+0.5*stepX,z-stepZ);
    double HminusPrev = eq->H(xPrevShifted-0.5*stepX,z-stepZ);
    double gvalPrev = eq->G(xPrevShifted,z-stepZ);
    double fvalPrev = eq->F(xPrevShifted,z-stepZ);

    // Fill right hand side
    cdouble left, center, right;
    if ( bTr != NULL )
    {
      bTr->threePointStencil( ix, iz, left, center, right);
    }
    else
    {
      //left = ix > 0 ? (*prevSolution)(ix-1):guide->transverseBC(z-stepZ, WaveGuideFDSimulation::Boundary_t::BOTTOM);
      left = ix > 0 ? (*prevSolution)(ix-1):bc->fixedField(x,z);
      center = (*prevSolution)(ix);
      //right = ix < Nx-1 ? (*prevSolution)(ix+1):guide->transverseBC(z-stepZ, WaveGuideFDSimulation::Boundary_t::TOP);
      right = ix < Nx-1 ? (*prevSolution)(ix+1):bc->fixedField(x,z);
    }

    if ( ix > 0 )
    {
      //rhs[ix] = (*solution)(ix-1,iz-1)*Hminus*gval;
      rhs[ix] = left*Hminus*gval;
    }
    else
    {
      //rhs[ix] = guide->transverseBC(z, WaveGuideFDSimulation::Boundary_t::BOTTOM)*Hminus*gval +\
      guide->transverseBC(z-stepZ, WaveGuideFDSimulation::Boundary_t::BOTTOM)*HminusPrev*gvalPrev; // Make sure that it is not a random value
      rhs[ix] = bc->fixedField(x,z)*Hminus*gval + left*HminusPrev*gvalPrev;

      diag[ix] -= 0.25*IMAG_UNIT*rho*bc->neighbourCoupling( *this, x, z );
    }

    if ( ix < Nx-1 )
    {
      //rhs[ix] += (*solution)(ix+1,iz-1)*HplussPrev*gvalPrev;
      rhs[ix] += right*HplussPrev*gvalPrev;
    }
    else
    {
      //rhs[ix] += ( guide->transverseBC(z, WaveGuideFDSimulation::Boundary_t::TOP)*Hpluss*gval + \
      guide->transverseBC(z-stepZ, WaveGuideFDSimulation::Boundary_t::TOP)*HplussPrev*gvalPrev );
      rhs[ix] += ( bc->fixedField(x,z)*Hpluss*gval + \
      right*HplussPrev*gvalPrev );

      diag[ix] -= 0.25*IMAG_UNIT*rho*bc->neighbourCoupling( *this, x, z );
    }

    rhs[ix] *=  (0.25*IMAG_UNIT*rho);
    //rhs[ix] -= 0.25*(*solution)(ix,iz-1)*IMAG_UNIT*rho*(HplussPrev+HminusPrev)*gval;
    rhs[ix] -= 0.25*center*IMAG_UNIT*rho*(HplussPrev+HminusPrev)*gval;
    //rhs[ix] += (1.0*fvalPrev - 0.5*(betaPrev*r + IMAG_UNIT*deltaPrev*r) )*(*solution)(ix,iz-1);
    rhs[ix] += (1.0*fvalPrev - 0.5*(betaPrev*r + IMAG_UNIT*deltaPrev*r) )*center;
  }

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
