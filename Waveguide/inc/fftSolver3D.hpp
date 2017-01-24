#ifndef FFT_SOLVER_3D_H
#define FFT_SOLVER_3D_H
#include "solver3D.hpp"
#include <complex>
#include <fftw3.h>

typedef std::complex<double> cdouble;
class FFTSolver3D: public Solver3D
{
public:
  FFTSolver3D(): Solver3D("FFTSolver3D"){};
  virtual ~FFTSolver3D();

  /** Propgate one step */
  void solveStep( unsigned int step );
private:
  /** The convolution kernel */
  cdouble kernel( double kx, double ky ) const;

  /** Return the spatial frequency corresponding to indx */
  double spatialFreq( unsigned int indx, unsigned int size ) const;

  fftw_complex *curr{NULL};
  fftw_complex *prev{NULL};
  fftw_plan ftforw;
  fftw_plan ftback;
  bool planInitialized{false};

  /** Diffraction step */
  void propagate();

  /** Refraction step */
  void refraction( unsigned int step );
};
#endif
