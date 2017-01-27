#include "postProcessMod.hpp"
#include "solver.hpp"
#include <iostream>

using namespace std;

void post::Intensity::result( const Solver &solver, arma::mat &res )
{
  res = arma::abs( solver.getSolution() );
}

void post::Intensity::result( const Solver &solver, arma::cube &res )
{
  res = arma::abs( solver.getSolution3D() );
}

void post::Phase::result( const Solver &solver, arma::mat &res )
{
  res = arma::arg( solver.getSolution() );
}

void post::Phase::result( const Solver &solver, arma::cube &res )
{
  res = arma::arg( solver.getSolution3D() );
}

void post::ExitField::result( const Solver &solver, arma::vec &res )
{
  res = arma::real( solver.getSolution().col(solver.getSolution().n_cols-1) );
}

void post::ExitField::result( const Solver &solver, arma::mat&res )
{
  res = arma::real( solver.getSolution3D().slice(solver.getSolution3D().n_slices-1) );
}

void post::ExitIntensity::result( const Solver &solver, arma::vec &res )
{
  res = arma::pow( arma::abs( solver.getSolution().col(solver.getSolution().n_cols-1) ), 2 );
}

void post::ExitIntensity::result( const Solver &solver, arma::mat &res )
{
  res = arma::pow( arma::abs( solver.getSolution3D().slice(solver.getSolution3D().n_slices-1) ), 2 );
}

void post::ExitPhase::result( const Solver &solver, arma::vec &res )
{
  res = arma::arg( solver.getSolution().col(solver.getSolution().n_cols-1) );
}

void post::ExitPhase::result( const Solver &solver, arma::mat &res )
{
  res = arma::arg( solver.getSolution3D().slice(solver.getSolution3D().n_slices-1) );
}
