#ifndef BEND_OPTIMIZER_H
#define BEND_OPTIMIZER_H
#include "multipleCurvedWG.hpp"
#include <map>
#include <string>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

class BendOptimizer
{
public:
  BendOptimizer(){};
  ~BendOptimizer();

  /** Initialize the simulation */
  void init( const std::map<std::string, double> &parameters );
private:
  MultipleCurvedWG wgs;
  Json::Value geometry;
  const std::map<std::string, double> *params;
  const gsl_multimin_fminimizer_type *T{ gsl_multimin_fminimizer_nmsimplex2 };
  gsl_multimin_fminimizer *minimizer{NULL};
  gsl_vector *variables{NULL};
  gsl_multimin_function optFunction;
  static double targetFunction( const gsl_vector *variables, void *par );

  /** Populates a JSON object based on the values in the GSL vector */
  void populateJSON( const gsl_vector *var );

  /** Populates a GSL vector based on the values in the JSON object */
  void populateGSLVector();
};
#endif
