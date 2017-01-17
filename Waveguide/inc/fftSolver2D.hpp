#ifndef FFT_SOLVER_2D_H
#define FFT_SOLVER_2D_H
#include "solver2D.hpp"
#include <complex>

typedef std::complex<double> cdouble;

class FFTSolver2D: public Solver2D
{
public:
  FFTSolver2D(): Solver2D("FFTSolver2D"){};
  virtual ~FFTSolver2D();

  void solveStep( unsigned int step ) override;
protected:
  unsigned int initialLength{0};

  /** Void init */
  void init();

  bool callInit{true};

  arma::cx_vec *filterSignal{new arma::cx_vec()};

  /** Kernel function */
  cdouble kernel( double kx ) const;

  /** Propagate an FFT step */
  void propagate();

  /** Perform the refraction step */
  void refraction( unsigned int step );

  /** Return the spatial frequency corresponding to indx */
  double spatialFreq( unsigned int indx, unsigned int size ) const;

  /** Return the part of the signal that should be included in the final solution */
  virtual arma::cx_vec& signalToFilter() const override;
};
#endif
