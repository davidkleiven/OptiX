#include "fftSolver3D.hpp"
#include <cmath>
#include <visa/visa.hpp>
#include "paraxialSimulation.hpp"
using namespace std;

const double PI = acos(-1.0);
typedef visa::Colormaps::Colormap_t cmap_t;
FFTSolver3D::~FFTSolver3D()
{
  if ( planInitialized )
  {
    fftw_destroy_plan( ftforw );
    fftw_destroy_plan( ftback );
  }
}

cdouble FFTSolver3D::kernel( double kx, double ky ) const
{
  double dz = guide->longitudinalDiscretization().step;
  double k = guide->getWavenumber();
  cdouble A(0.0, 0.5*dz/k);
  return exp(-A*(kx*kx + ky*ky) );
}

double FFTSolver3D::spatialFreq( unsigned int indx, unsigned int size ) const
{
  if ( indx > size/2 )
  {
    indx = size-indx;
  }
  double stepX = guide->transverseDiscretization().step;
  return 2.0*PI*indx/(stepX*size);
}

void FFTSolver3D::solveStep( unsigned int step )
{
  if ( !planInitialized )
  {
    assert( prevSolution->is_square() );
    assert( currentSolution->is_square() );

    prev = reinterpret_cast<fftw_complex*>( prevSolution->memptr() );
    curr = reinterpret_cast<fftw_complex*>( currentSolution->memptr() );
    ftforw = fftw_plan_dft_2d( prevSolution->n_rows, prevSolution->n_cols, prev, curr, FFTW_FORWARD, FFTW_ESTIMATE );
    ftback = fftw_plan_dft_2d( prevSolution->n_rows, prevSolution->n_cols, curr, prev, FFTW_BACKWARD, FFTW_ESTIMATE );
    planInitialized = true;
  }
  propagate();
  refraction( step );
}

void FFTSolver3D::propagate()
{
  // NOTE: FFTW assumes row-major ordering, while Armadillo uses column major
  fftw_execute( ftforw ); // prev --> current

  if ( visFourierSpace )
  {
    arma::mat fourierIntensity = arma::abs(*currentSolution );
    plots.get("FourierIntensity").fillVertexArray( fourierIntensity );
    plots.show();
  }

  for ( unsigned int i=0;i<currentSolution->n_cols;i++ )
  {
    double kx = spatialFreq( i, currentSolution->n_cols );
    for ( unsigned int j=0;j<currentSolution->n_rows;j++ )
    {
      double ky = spatialFreq( j, currentSolution->n_rows );
      (*currentSolution)(j,i) *= kernel( kx, ky );
    }
  }

  fftw_execute( ftback ); // current --> prev
}

void FFTSolver3D::refraction( unsigned int step )
{
  assert ( step > 0 );
  double stepZ = guide->longitudinalDiscretization().step;
  double wavenumber = guide->getWavenumber();
  double z1 = guide->getZ( step );
  double z0 = guide->getZ( step-1 );
  cdouble im(0.0,1.0);
  const double ZERO = 1E-10;
  for ( unsigned int i=0;i<prevSolution->n_cols; i++ )
  {
    double x = guide->getX(i);
    for ( unsigned int j=0;j<prevSolution->n_rows;j++ )
    {
      double y = guide->getX(j);
      double delta, beta, deltaPrev, betaPrev;
      guide->getXrayMatProp( x, y, z1, delta, beta );
      guide->getXrayMatProp( x, y, z0, deltaPrev, betaPrev );
      if (( abs(delta-deltaPrev) < ZERO ) || ( abs(beta-betaPrev) < ZERO ))
      {
        // Wave has crossed a border
        refractionIntegral( x, y, z0, z1, delta, beta );
      }

      // FFTW3: Divide by length to normalize
      double normalization = prevSolution->n_rows*prevSolution->n_cols;
      //normalization = 1.0;
      (*currentSolution)(j,i) = (*prevSolution)(j,i)*exp( -wavenumber*(beta+im*delta)*stepZ )/normalization;
    }
  }


  if ( visRealSpace )
  {
    arma::mat values = arma::abs( *currentSolution );
    plots.get("Intensity").fillVertexArray( values );
    values = -arma::arg( *currentSolution );
    plots.get("Phase").fillVertexArray( values );
    plots.show();
  }
}

void FFTSolver3D::visualizeRealSpace()
{
  visRealSpace = true;
  plots.addPlot("Intensity");
  plots.addPlot("Phase");
  plots.get("Intensity").setCmap( cmap_t::NIPY_SPECTRAL );
  plots.get("Phase").setCmap( cmap_t::NIPY_SPECTRAL );
}

void FFTSolver3D::visualizeFourierSpace()
{
  visFourierSpace = true;
  plots.addPlot("FourierIntensity");
}

void FFTSolver3D::setIntensityMinMax( double min, double max )
{
  plots.get("Intensity").setColorLim( min, max );
}

void FFTSolver3D::setPhaseMinMax( double min, double max )
{
  plots.get("Phase").setColorLim( min, max );
}

void FFTSolver3D::refractionIntegral( double x, double y, double z1 , double z2, double &delta, double &beta )
{
  double deltaTemp = 0.0;
  double betaTemp = 0.0;
  guide->getXrayMatProp( x, y, z1, delta, beta );
  guide->getXrayMatProp( x, y, z2, deltaTemp, betaTemp );
  delta += deltaTemp;
  beta += betaTemp;
  double z = z1;
  double dz = (z2-z1)/nStepsInRefrIntegral;
  for ( unsigned int i=1;i<nStepsInRefrIntegral-1;i++ )
  {
    z = z1 + i*dz;
    guide->getXrayMatProp( x, y, z, deltaTemp, betaTemp );
    delta += 2.0*deltaTemp;
    beta += 2.0*betaTemp;
  }
  // NOTE: Multiplication with dz is left out as this is done in the refraction function
  delta /= (2.0*nStepsInRefrIntegral);
  beta /= (2.0*nStepsInRefrIntegral);
}
