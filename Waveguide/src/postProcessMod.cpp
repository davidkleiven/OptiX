#include "postProcessIntensity.hpp"

void post::Intensity::result( const Solver2D &solver, arma::mat &res )
{
  res = arma::abs( solver.getSolution() );
}

void post::Phase::result( const Solver2D &solver, arma::mat &res )
{
  res = arma::arg( solver.getSolution() );
}
