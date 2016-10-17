#include "straightWG2D.hpp"
#include "controlFile.hpp"
#include "solver2D.hpp"

using namespace std;
bool StraightWG2D::isInsideGuide( double x, double z ) const
{
  return ( x > 0.0 ) && ( x < width );
}

void StraightWG2D::fillInfo( Json::Value &obj ) const
{
    obj["Width"] = width;
}

double StraightWG2D::waveGuideStartX( double z ) const
{
  return 0.0;
}

double StraightWG2D::waveGuideEndX( double z ) const
{
  return width;
}

void StraightWG2D::init( const ControlFile &ctl )
{
  WaveGuideFDSimulation::init(ctl);
  width = ctl.get()["waveguide"]["Width"].asDouble();
}

void StraightWG2D::extractField( double wcrd, vector<cdouble> &res ) const
{
  double z0 = 0.0;
  unsigned int ix, dummy;
  res.clear(); // Make sure that the vector is empty
  closestIndex( wcrd, z0, ix, dummy );
  for ( unsigned int iz=0; iz<nodeNumberLongitudinal(); iz++ )
  {
    res.push_back( solver->getSolution()(ix, iz) );
  }
}
