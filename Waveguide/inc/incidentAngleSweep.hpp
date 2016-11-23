#ifndef INCIDENT_ANGLE_SWEEP_H
#define INCIDENT_ANGLE_SWEEP_H
#include "straightWG2D.hpp"
#include "planeWave.hpp"
#include "cladding.hpp"
#include "paraxialEquation.hpp"
#include "crankNicholson.hpp"
#include "visualizer.hpp"
#include <armadillo>
#include <set>
#include <string>

/**
* Class for running a sweep over incident angle and saving the farfields
*/
class IncidentAngleSweep
{
public:
  IncidentAngleSweep();
  void setWidth( double w );
  void setWavelength( double wl );
  void setLength( double L );
  void setIncAngles( double min, double max, unsigned int Nangles );
  void save( const std::string &fname ) const;
  void solve();
  void setTransverseDisc( double xmin, double xmax, unsigned int Nx );
  void setLongitudinalDisc( double zmin, double zmax, unsigned int Nz );
  void setCladdingSilicon( double energyInEv );
  void setAlcoholInside( double energyInEv );
  void setEthylenGlycolInside( double energyInEv );
  void saveIndx( unsigned int indx );
  void setFFTSignalLength( unsigned int length){fftSignalLength=length;};
  void setVisualizer( Visualizer &newvis ){ vis = &newvis; };
  void savePic( const char *dir );
private:
  double getAngle( unsigned int indx ) const;
  StraightWG2D wg;
  PlaneWave pw;
  Cladding cladding;
  Cladding inside;
  ParaxialEquation eq;
  CrankNicholson solver;
  arma::mat farField;
  double thetaMin{0.0};
  double thetaMax{0.0};
  unsigned int nTheta{1};
  std::set<unsigned int> indxToSave;
  unsigned int uid{0};
  unsigned int fftSignalLength{65536};
  Visualizer *vis{NULL};
  std::string picDir;
  bool saveImages{false};
};
#endif
