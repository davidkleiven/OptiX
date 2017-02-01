#include "solver.hpp"
#include <stdexcept>

using namespace std;
const arma::cx_cube& Solver::getSolution3D() const
{
  throw ( runtime_error("The 3D version of getSolution() is not implemented!") );
}

const arma::cx_mat& Solver::getSolution() const
{
  throw ( runtime_error("The 2D version of getSolution() is not implemented!") );
}

const arma::cx_vec& Solver::getLastSolution() const
{
  throw ( runtime_error("The 2D version of getLastSolution() is not implemented!") );
}

const arma::cx_mat& Solver::getLastSolution3D() const
{
  throw( runtime_error("The 3D version of getLastSolution3D() is not implemented!") );
}

arma::cx_mat& Solver::getLastSolution3D()
{
  throw( runtime_error("The 3D version of non-const getLastSolution3D() is not implemented!") );
}
