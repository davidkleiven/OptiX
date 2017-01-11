#ifndef FRESNEL_PROPAGATOR_H
#define FRESNEL_PROPAGATOR_H
#include <armadillo>
#include <complex>
#include <stdexcept>
#include "waveGuideFDSimulation.hpp"
typedef std::complex<double> cdouble;

/** Class for evolving an sourece according to the Fresnel convolution integral */
class FresnelPropagator
{
public:
  FresnelPropagator(){};
  FresnelPropagator( const FresnelPropagator &other ) = delete;
  FresnelPropagator& operator =( const FresnelPropagator other ) = delete;
  ~FresnelPropagator();

  /** Propagates the field a Nsteps */
  void propagate( unsigned int Nsteps );

  /** Set the stepsize in the forward direction in nano meters */
  void setStepsize( double step ){ dz = step; };

  /** Set the wavelength in nano meter */
  void setWavelength( double lambda );

  /** Set the initial field. The template parameter has to implement operator( double x ) that evaluates the field */
  template <class T>
  void setInitialConditions( const T &initfield );

  /** Gaussian kernel in the Fourier domain */
  cdouble kernel( double kx ) const;

  /** Sets the transverse discretization */
  void setTransverseDiscretization( double xmin, double xmax, unsigned int Nsteps );

  /** Saves the field to a HDF5 file */
  void save( const std::string &fname ) const;

  /** Set the pad factor. Size of the padded signal is padFactor*2^{floor(log2(original signal size))} */
  void setPadFactor( unsigned int factor ){ padFactor = factor; };

  /** Get the resulting intensity */
  const arma::mat& getIntensity() const { return *intensity; };
private:
  void step();

  /** Returns the spatial frequency corresponding to the array index */
  double spatialFreq( unsigned int indx, unsigned int size ) const;

  /** Get the x-coordinate corresponding to the array index ix */
  double getX( unsigned int ix ) const;

  /** Shift the FFT to have the zero frequency in the center */
  static void fftshift( arma::cx_vec &vec );

  /** Perform zero padding before doing FFT */
  void padSignal( const arma::cx_vec &vec, arma::cx_vec &padded ) const;

  /** Unpad the signal by extracting the central part of the padded signal */
  void unpadSignal( const arma::cx_vec &padded, arma::cx_vec &vec ) const;

  arma::mat *intensity{NULL};
  arma::cx_vec prev;
  double dz{0.0};
  double wavenumber{0.0};
  Disctretization *xDisc{NULL};
  unsigned int current{0};
  bool initConditionsSet{false};
  unsigned int padFactor{4};
};

// Implement hte template set function
template <class T>
void FresnelPropagator::setInitialConditions( const T &initfield )
{
  if ( xDisc == NULL )
  {
    throw (std::runtime_error("No transverse discretization is set!"));
  }

  for ( unsigned int ix=0;ix<prev.n_elem;ix++ )
  {
    double x = getX(ix);
    prev(ix) = initfield(x);
  }
  initConditionsSet = true;
};
#endif
