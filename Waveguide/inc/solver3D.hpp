#ifndef SOLVER_3D_H
#define SOLVER_3D_H
#include <armadillo>
#include "solver.hpp"

class ParaxialSimulation;

class Solver3D: public Solver
{
public:
  Solver3D( const char* name ):Solver(name, Dimension_t::THREE_D){};

  /** Set the simulation routine */
  void setSimulator( ParaxialSimulation &sim ) override;

  /** Set the initial conditions */
  void setInitialConditions( const arma::cx_mat &values ) override;

  /** Propagate one step */
  virtual void step() override;

  /** Solve the entire system */
  virtual void solve() override;
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

  /** Solver specific function to propagate one step */
  virtual void solveStep( unsigned int step ) = 0;
};
#endif
