#include "transmittivity.hpp"
#include "curvedWaveGuide2D.hpp"
#include "solver2D.hpp"
#include <cassert>
#include <cmath>

using namespace std;
post::Transmittivity::Transmittivity(): transmission(new vector<double>()){};

post::Transmittivity::~Transmittivity()
{
  delete transmission;
}

double post::Transmittivity::trapezoidalIntegrateIntensityZ( unsigned int ixStart, unsigned int ixEnd ) const
{
  double integral = pow( abs(guide->getSolver().getLastSolution()(ixStart) ),2) +\
    pow( abs(guide->getSolver().getLastSolution()(ixEnd) ),2);
  for ( unsigned int ix=ixStart+1; ix <= ixEnd-1; ix ++ )
  {
    double x = guide->getX(ix);
    integral += 2.0*pow( abs(guide->getSolver().getLastSolution()(ix)), 2 );
  }
  double dx = ( guide->getX( ixEnd) - guide->getX( ixStart ) )/static_cast<double>( ixEnd - ixStart );
  return integral*dx*0.5;
}

void post::Transmittivity::compute( double z )
{
  assert( guide != NULL );

  unsigned int wgStart = 0;
  unsigned int wgEnd = 0;
  unsigned int zIndx;
  
  // Integrate across waveguide
  if ( computeIntensityAtZero )
  {
    guide->closestIndex( 0.0, 0.0, wgStart, zIndx );
    guide->closestIndex( guide->getWidth(), 0.0, wgEnd, zIndx );
    intensityAtZero = trapezoidalIntegrateIntensityZ( wgStart, wgEnd );
    computeIntensityAtZero = false;
  }


  double xWgStart, xWgEnd;
  if ( guide->getBorderTracker() == NULL )
  {
    xWgStart = guide->waveGuideStartX( z );
    xWgEnd = guide->waveGuideEndX( z );

    // Assertions for debugging (allow first and last point to be outside )
    assert( guide->isInsideGuide( xWgStart+1, z ) );
    assert( guide->isInsideGuide( xWgEnd-1, z ) );
  }
  else
  {
    // Reset wg start and wgEnd as these should be constant throughout the waveguide if the border tracker is used
    xWgStart = guide->waveGuideStartX( 0.0 );
    xWgEnd = guide->waveGuideEndX( 0.0 );
  }

  guide->closestIndex( xWgStart, z, wgStart, zIndx );
  guide->closestIndex( xWgEnd, z, wgEnd, zIndx );

  double intensity = trapezoidalIntegrateIntensityZ( wgStart, wgEnd );
  transmission->push_back( intensity/intensityAtZero );
}

void post::Transmittivity::linkWaveguide( const CurvedWaveGuideFD &wg )
{
  guide = &wg;
}
