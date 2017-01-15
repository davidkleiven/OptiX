#include "bendOptimizer.hpp"
#include <cassert>

using namespace std;

BendOptimizer::~BendOptimizer()
{
  if ( variables != NULL ) gsl_vector_free( variables );
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
  variables = gsl_vector_alloc( 2*params->at("nWaveguides")-1 );
  populateGSLVector();

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
    double radius = gsl_vector_get( vec, i );
    radius = radius > 1E-5 ? radius:1E-5;

    geometry["waveguides"][i]["radius"] = radius;
    if ( i < nGuides-1 )
    {
      double angle =  gsl_vector_get( vec, i+nGuides );
      angle = angle > 0.0 ? angle:0.0;
      geometry["waveguides"][i]["angle"] = angle;
      totAngle += angle;
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

double BendOptimizer::targetFunction( const gsl_vector *vec, void *par )
{
  BendOptimizer* self = static_cast<BendOptimizer*>( par );
  self->populateJSON( vec );
  self->wgs.loadWaveguides( self->geometry );
  self->wgs.init( *self->params );
  self->wgs.solve();

  // Maxmimize transmittivity <=> minimize -transmittivity
  return -self->wgs.getTransmittivity()( self->wgs.getTransmittivity().n_elem-1 );
}

void BendOptimizer::optimize()
{
  gsl_multimin_function optFunction;
  optFunction.f = targetFunction;
  optFunction.n = variables->size;
  optFunction.params = this;
  gsl_multimin_fminimizer *minimizer = gsl_multimin_fminimizer_alloc( T, variables->size );
  gsl_vector *stepsize = gsl_vector_alloc( variables->size );
  for ( unsigned int i=0;i<variables->size; i++ )
  {
    gsl_vector_set( stepsize, i, 1.0 );
  }

  gsl_multimin_fminimizer_set( minimizer, &optFunction, variables, stepsize );

  for ( unsigned int iter=0; iter<params->at("maxIter");iter++ )
  {
    int status = gsl_multimin_fminimizer_iterate( minimizer );

    if ( status == GSL_SUCCESS )
    {
      clog << "Converged to minimum!\n";
      return;
    }
    transmittivity.push_back( -gsl_multimin_fminimizer_minimum( minimizer ) );
  }

  gsl_multimin_fminimizer_free( minimizer );
  gsl_vector_free( stepsize );
  clog << "The maximum number of iterations was reached!\n";
}
