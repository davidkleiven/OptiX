#ifndef FFT_SOLVER_3D_H
#define FFT_SOLVER_3D_H
#include "solver3D.hpp"
#include <complex>
#include <fftw3.h>
#include <visa/visa.hpp>

typedef std::complex<double> cdouble;
class FFTSolver3D: public Solver3D
{
public:
  FFTSolver3D(): Solver3D("FFTSolver3D"){};
  virtual ~FFTSolver3D();

  /** Propgate one step */
  void solveStep( unsigned int step );

  /** A call to this before solving will visualize the real space intensity on each iteration */
  void visualizeRealSpace();

  /** A call to this before solving will visualize the fourier intensity on each iteration */
  void visualizeFourierSpace();

  /** Set the maximum and minimum value in the intensity visualization */
  void setIntensityMinMax( double min, double max );

  /** Set the maximum and minimum value in the phase visualization */
  void setPhaseMinMax( double min, double max );

  /** Resets the solver */
  virtual void reset() override;
private:
  /** The convolution kernel */
  cdouble kernel( double kx, double ky ) const;

  /** Return the spatial frequency corresponding to indx in x-direction */
  double spatialFreqX( unsigned int indx, unsigned int size ) const;

  /** Return the spatial frequency corresponding to indx in the y-direction */
  double spatialFreqY( unsigned int indx, unsigned int size ) const;

  fftw_complex *curr{NULL};
  fftw_complex *prev{NULL};
  fftw_plan ftforw;
  fftw_plan ftback;
  visa::WindowHandler plots;
  bool planInitialized{false};

  /** Diffraction step */
  void propagate();

  /** Refraction step */
  void refraction( unsigned int step );

  /** Computes the refraction integral when a border has been crossed */
  void refractionIntegral( double x, double y, double z1, double z2, double &delta, double &beta );

  /** Evaluates the refractive index for the purpose of overlay */
  void evaluateRefractiveIndex( arma::mat &refr, double z ) const;

  unsigned int nStepsInRefrIntegral{200};

  bool visRealSpace{false};
  bool visFourierSpace{false};
};

class FFT3DSolverDebug: public FFTSolver3D
{
public:
  /** Only for debugging. The non-const version should not be used */
  virtual arma::cx_mat& getLastSolution3D() override{ return *prevSolution; };
};
#endif
