#include "fresnelPropagator.hpp"
#include <cmath>
#include <H5Cpp.h>
#include <hdf5_hl.h>
//#define DEBUG_FFT

using namespace std;

const double PI = acos(-1.0);

cdouble FresnelPropagator::kernel( double kx ) const
{
  cdouble A(0.0, 0.5/wavenumber);
  return exp( -A*dz*kx*kx);
}

FresnelPropagator::~FresnelPropagator()
{
  if ( xDisc != NULL ) delete xDisc;
  if ( intensity != NULL ) delete intensity;
}

void FresnelPropagator::setWavelength( double lambda )
{
  wavenumber = 2.0*PI/lambda;
}

void FresnelPropagator::setTransverseDiscretization( double xmin, double xmax, unsigned int nSteps )
{
  if ( xDisc == NULL )
  {
    xDisc = new Disctretization();
  }

  unsigned int lg2 = log2(nSteps);
  if ( pow(2, lg2) != nSteps )
  {
    clog << "Warning: the algorithm is most efficient if the number of steps is an integer power of 2\n";
  }
  xDisc->min = xmin;
  xDisc->max = xmax;
  xDisc->step = (xmax-xmin)/nSteps;
  prev.set_size(nSteps);
}

double FresnelPropagator::spatialFreq( unsigned int indx ) const
{
  unsigned int N = intensity->n_rows;
  if ( indx > N/2 )
  {
    indx = N-indx;
  }
  return 2.0*PI*indx/(xDisc->step*N);// - PI/(xDisc->step*N);
}

void FresnelPropagator::step()
{
  // Fill in the amplitude
  for ( unsigned int i=0;i<intensity->n_rows;i++ )
  {
    (*intensity)(i,current) = abs( prev(i) );
  }

  // Fourier transform previous step
  arma::cx_vec fftPrev = arma::fft( prev );

  #ifdef DEBUG_FFT
    if ( current == 0 )
    {
      arma::vec vec = arma::abs(fftPrev);
      string fname("data/fftDebug.h5");
      vec.save(fname.c_str(), arma::hdf5_binary);
      clog << "FFT debug written to " << fname << endl;
    }
  #endif

  for ( unsigned int i=0;i<intensity->n_rows;i++ )
  {
    double kx = spatialFreq(i);
    fftPrev(i) *= kernel(kx);
  }

  // Inverse fourier transform
  prev = arma::ifft(fftPrev);
  //fftshift(prev);
  current++;
}

void FresnelPropagator::save( const string &fname ) const
{
  int uid = rand()%10000000;
  if ( intensity == NULL )
  {
    throw (runtime_error("No solution has been computed!"));
  }

  hsize_t dim[2] = {intensity->n_cols, intensity->n_rows};
  stringstream ss;
  ss << fname << uid << ".h5";
  hid_t file_id = H5Fcreate( ss.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  H5LTmake_dataset( file_id, "intensity", 2, dim, H5T_NATIVE_DOUBLE, intensity->memptr());

  // Set useful attributes
  H5LTset_attribute_double( file_id, "intensity", "wavenumber", &wavenumber, 1);
  H5LTset_attribute_double( file_id, "intensity", "xmin", &xDisc->min, 1);
  H5LTset_attribute_double( file_id, "intensity", "xmax", &xDisc->max, 1);
  H5LTset_attribute_double( file_id, "intensity", "xstep", &xDisc->step, 1);
  H5LTset_attribute_double( file_id, "intensity", "zstep", &dz, 1);
  double zmin = 0.0;
  H5LTset_attribute_double( file_id, "intensity", "zmin", &zmin, 1);
  double zmax = dz*intensity->n_cols;
  H5LTset_attribute_double( file_id, "intensity", "zmax", &zmax, 1);
  H5LTset_attribute_int( file_id, "intensity", "uid", &uid, 1);
  H5Fclose(file_id);

  clog << "Intensity written to " << ss.str() << endl;
}

double FresnelPropagator::getX( unsigned int ix ) const
{
  return xDisc->min + xDisc->step*ix;
}

void FresnelPropagator::propagate( unsigned int nSteps )
{
  if ( !initConditionsSet )
  {
    throw (runtime_error("Initial conditions are not set!"));
  }

  if ( intensity == NULL )
  {
    intensity = new arma::mat();
  }
  intensity->set_size( prev.n_elem, nSteps+1 );
  for ( unsigned int i=0;i<nSteps;i++ )
  {
    step();
  }
}

void FresnelPropagator::fftshift( arma::cx_vec &vec )
{
  for ( unsigned int i=0;i<vec.n_elem/2;i++ )
  {
    cdouble copy = vec(i);
    vec(i) = vec(i+vec.n_elem/2);
    vec(i+vec.n_elem/2) = copy;
  }
}
