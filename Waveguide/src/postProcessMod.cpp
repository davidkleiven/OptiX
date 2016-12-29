#include "postProcessMod.hpp"
#include "solver2D.hpp"

void post::Intensity::result( const Solver2D &solver, arma::mat &res )
{
  res = arma::abs( solver.getSolution() );
}

void post::Phase::result( const Solver2D &solver, arma::mat &res )
{
  res = arma::arg( solver.getSolution() );
}

void post::ExitField::result( const Solver2D &solver, arma::vec &res )
{
  res = arma::real( solver.getSolution().col(solver.getSolution().n_cols-1) );
}

void post::ExitIntensity::result( const Solver2D &solver, arma::vec &res )
{
  res = arma::abs( solver.getSolution().col(solver.getSolution().n_cols-1) );
}
