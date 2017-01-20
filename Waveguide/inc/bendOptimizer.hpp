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

  /** Saves a JSON file with the optimal parameters */
  void save( const char* fname ) const;
private:
  MultipleCurvedWG wgs;
  Json::Value geometry;
  const std::map<std::string, double> *params;
  const gsl_multimin_fminimizer_type *T{ gsl_multimin_fminimizer_nmsimplex2 };
  gsl_vector *variables{NULL};
  static double targetFunction( const gsl_vector *variables, void *par );
  std::vector<double> arclengths;
  double prevRunResult{1E-16};

  /** Populates a JSON object based on the values in the GSL vector */
  void populateJSON( const gsl_vector *var );

  /** Populates a GSL vector based on the values in the JSON object */
  void populateGSLVector();
  std::vector<double> transmittivity; // Track the transmittivity to study the sensitivity

  /** Computes the accumulated relative change in the arclength */
  double accumulatedRelativeArcLengthChange() const;

  /** Computes the arclength */
  static double arcLength( double r, double thetaDeg );

  /** Generates a random number in the interval [0,1] */
  static double uniform();

  /** Adjusts the poposed configuration vectors to the constraints */
  void adjustToConstraints( gsl_vector *vec ) const;
};
#endif
