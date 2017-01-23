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
