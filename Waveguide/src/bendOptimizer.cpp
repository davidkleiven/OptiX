#include "bendOptimizer.hpp"
#include <cassert>

using namespace std;

BendOptimizer::~BendOptimizer()
{

}
void BendOptimizer::init( const map<string,double> &parameters )
{
  params = &parameters;
  geometry["waveguides"] = Json::Value( Json::arrayValue );
  for ( unsigned int i=0;i<params->at("nWaveguides");i++ )
  {
    Json::Value val;
    val["radius"] = 40.0;
    val["curvature"] = "concave";
    val["angle"] = params->at("totalDeflectionAngle")/params->at("nWaveguides");
    geometry["waveguides"].append( val );
  }

  wgs.loadWaveguides( geometry );
  wgs.init( *params );
}

void BendOptimizer::populateJSON( const gsl_vector *vec )
{
  unsigned int nGuides = geometry["waveguides"].size();
  assert( vec->size == 2*nGuides-1 );

  double totAngle = 0.0;
  for ( unsigned int i=0;i<nGuides;i++ )
  {
    geometry["waveguides"][i]["radius"] = gsl_vector_get( vec, i );
    if ( i < nGuides-1 )
    {
      geometry["waveguides"][i]["angle"] = gsl_vector_get( vec, i+nGuides );
      totAngle += geometry["waveguides"][i]["angle"].asDouble();
    }
    else
    {
      geometry["waveguides"][i]["angle"] = params->at("totalDeflectionAngle") - totAngle;
      assert( params->at("totalDeflectionAngle") - totAngle >= 0.0 );
    }
  }
}

void BendOptimizer::populateGSLVector()
{
  unsigned int nGuides = geometry["waveguides"].size();
  assert( variables->size == 2*nGuides-1 );
  for ( unsigned int i=0;i<nGuides;i++ )
  {
    gsl_vector_set( variables, i, geometry["waveguides"][i]["radius"].asDouble() );

    if ( i < nGuides-1 )
    {
      gsl_vector_set( variables, i+nGuides, geometry["waveguides"][i]["angle"].asDouble() );
    }
  }
}
