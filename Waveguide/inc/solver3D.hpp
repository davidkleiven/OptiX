#ifndef SOLVER_3D_H
#define SOLVER_3D_H
#include <armadillo>
#include "solver.hpp"

class ParaxialSimulation;

class Solver3D: public Solver
{
public:
  Solver3D( const char* name ):Solver(name, Dimension_t::THREE_D){};
  void setSimulator( ParaxialSimulation &sim ) override;

protected:
  arma::cx_cube *solution{NULL};
  arma::cx_mat *currentSolution{NULL};
  arma::cx_mat *prevSolution{NULL};

  // Some parameters
  unsigned int Nx, Ny, Nz;

  /** Filter the transverse signal */
  void filterTransverse( arma::cx_mat &mat );

  /** Copy the current solution to the previous array */
  void copyCurrentSolution( unsigned int step );
};
#endif
