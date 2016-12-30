#ifndef CRANC_NICHOLSON_H
#define CRANC_NICHOLSON_H
#include "solver2D.hpp"
#include "thomasAlgorithm.hpp"
#include <complex>

template<class T>
struct CSRBuild
{
  T *elements{NULL};
  T *row{NULL};
  T *col{NULL};
};

/** Class for handling the Crank-Nicholson discretization scheme */
class CrankNicholson: public Solver2D
{
public:
  CrankNicholson():Solver2D("CrankNicholson"){};
  ~CrankNicholson();

  /** Solves single step */
  virtual void solveStep( unsigned int iz ) override final;
protected:
  ThomasAlgorithm matrixSolver;

  /** Performs one iteration */
  void solveCurrent( unsigned int iz );
};

#endif
