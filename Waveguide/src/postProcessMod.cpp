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
  res = arma::real( solver.getLastSolution() );
}

void post::ExitField::result( const Solver &solver, arma::mat&res )
{
  res = arma::real( solver.getLastSolution3D() );
}

void post::ExitIntensity::result( const Solver &solver, arma::vec &res )
{
  res = arma::pow( arma::abs( solver.getLastSolution() ), 2 );
}

void post::ExitIntensity::result( const Solver &solver, arma::mat &res )
{
  res = arma::pow( arma::abs( solver.getLastSolution3D() ), 2 );
}

void post::ExitPhase::result( const Solver &solver, arma::vec &res )
{
  res = arma::arg( solver.getLastSolution() );
}

void post::ExitPhase::result( const Solver &solver, arma::mat &res )
{
  res = arma::arg( solver.getLastSolution3D() );
}
