#ifndef FRESNEL_PROPAGATOR_H
#define FRESNEL_PROPAGATOR_H
#include <armadillo>
#include <complex>
#include <stdexcept>
typedef std::complex<double> cdouble;
class Disctretization;

class FresnelPropagator
{
public:
  FresnelPropagator(){};
  ~FresnelPropagator();
  void propagate( unsigned int Nsteps );
  void setStepsize( double step ){ dz = step; };
  void setWavelength( double lambda );
  template <class T>
  void setInitialConditions( const T &initfield );
  cdouble kernel( double kx ) const;
  void setTransverseDiscretization( double xmin, double xmax, unsigned int Nsteps );
  void save( const std::string &fname ) const;
private:
  void step();
  double spatialFreq( unsigned int indx ) const;
  double getX( unsigned int ix ) const;
  static void fftshift( arma::cx_vec &vec );
  arma::mat *intensity{NULL};
  arma::cx_vec prev;
  double dz{0.0};
  double wavenumber{0.0};
  Disctretization *xDisc{NULL};
  unsigned int current{0};
  bool initConditionsSet{false};
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
