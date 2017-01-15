#ifndef BEND_OPTIMIZER_H
#define BEND_OPTIMIZER_H
#include "multipleCurvedWG.hpp"
#include <map>
#include <string>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <vector>

class BendOptimizer
{
public:
  BendOptimizer(){};
  ~BendOptimizer();

  /** Initialize the simulation */
  void init( const std::map<std::string, double> &parameters );

  /** Run the optimization */
  void optimize();
private:
  MultipleCurvedWG wgs;
  Json::Value geometry;
  const std::map<std::string, double> *params;
  const gsl_multimin_fminimizer_type *T{ gsl_multimin_fminimizer_nmsimplex2 };
  gsl_vector *variables{NULL};
  static double targetFunction( const gsl_vector *variables, void *par );

  /** Populates a JSON object based on the values in the GSL vector */
  void populateJSON( const gsl_vector *var );

  /** Populates a GSL vector based on the values in the JSON object */
  void populateGSLVector();
  std::vector<double> transmittivity; // Track the transmittivity to study the sensitivity
};
#endif
