#ifndef INCIDENT_ANGLE_SWEEP_H
#define INCIDENT_ANGLE_SWEEP_H
#include "straightWG2D.hpp"
#include "planeWave.hpp"
#include "cladding.hpp"
#include "paraxialEquation.hpp"
#include "crankNicholson.hpp"
#include <armadillo>
#include <set>

/**
* Class for running a sweep over incident angle and saving the farfields
*/
class IncidentAngleSweep
{
public:
  IncidentAngleSweep(){};
  void setWidth( double w );
  void setWavelength( double wl );
  void setLength( double L );
  void setIncAngles( double min, double max, Nangles );
  void save( const std::string &fname ) const;
  void solve();
  void setTransverseDisc( double xmin, double xmax, unsigned int Nx );
  void setLongitudinalDisc( double zmin, double zmax, unsigned int Nz );
  void setCladdingSilicon();
  void saveIndx( unsigned int indx );
private:
  double getAngle( indx ) const;
  StraightWG wg;
  PlaneWave pw;
  Cladding cladding;
  ParaxialEquation eq;
  CrankNicholson solver;
  arma::mat farField;
  double theta_min{0.0};
  double theta_max{0.0};
  unsigned int nTheta{1};
  std::set<unsigned int> indxToSave;
};
#endif
