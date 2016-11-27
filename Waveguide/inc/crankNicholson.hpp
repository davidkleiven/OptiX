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

  /** Solves the system */
  void solve() override final;

  /** Set the required parameters from the waveguide object */
  void initValuesFromWaveGuide();
protected:
  unsigned int Nx{0};
  unsigned int Nz{0};
  double stepX{1.0}, stepZ{1.0};
  double xmin{0.0};
  double zmin{0.0};
  double wavenumber{1.0};
  ThomasAlgorithm matrixSolver;

  /** Performs one iteration */
  void solveCurrent( unsigned int iz );
};

#endif
