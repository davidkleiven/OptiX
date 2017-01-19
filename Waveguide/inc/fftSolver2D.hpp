#ifndef FFT_SOLVER_2D_H
#define FFT_SOLVER_2D_H
#include "solver2D.hpp"
#include <complex>

typedef std::complex<double> cdouble;

class FFTSolver2D: public Solver2D
{
public:
  FFTSolver2D(): Solver2D("FFTSolver2D"){};

  /** Propagate the beam one step */
  void solveStep( unsigned int step ) override;
protected:
  unsigned int initialLength{0};

  /** Kernel function */
  cdouble kernel( double kx ) const;

  /** Propagate an FFT step */
  void propagate();

  /** Perform the refraction step */
  void refraction( unsigned int step );

  /** Return the spatial frequency corresponding to indx */
  double spatialFreq( unsigned int indx, unsigned int size ) const;
};
#endif
