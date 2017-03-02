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

  /** Parameters used internally, but needs to be shared between member functions */
  double Hpluss, Hminus, gval, HplussPrev, HminusPrev, gvalPrev, rho;

  bool printBC{true};

  /** Performs one iteration */
  void solveCurrent( unsigned int iz );

  /** Apply boundary condition to the solution matrix */
  void applyBC( cdouble subdiag[], cdouble diag[], cdouble rhs[] );

  /** Apply transparent boundary conditions */
  void applyTBC( cdouble subdiag[], cdouble diag[], cdouble rhs[] );

  /** Apply TBC on the side of the matrix/right hand side */
  void applyTBCOneSide( cdouble subdiag[], cdouble diag[], cdouble rhs[], unsigned int outer, unsigned int inner );
};

#endif
