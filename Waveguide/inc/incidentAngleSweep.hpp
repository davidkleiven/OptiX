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

  /** Set the width of the waveguide in nano meters */
  void setWidth( double w );

  /** Set the vacuum wavelength in nano meters */
  void setWavelength( double wl );

  /** Set the length of the waveguide in milli meters */
  void setLength( double L );

  /** Set the incident angles to run */
  void setIncAngles( double min, double max, unsigned int Nangles );

  /** Save the results to a HDF5 file */
  virtual void save( const std::string &fname ) const;

  /** Run the simulation */
  void solve();

  /** Set the transverse discretization */
  void setTransverseDisc( double xmin, double xmax, unsigned int Nx );

  /** Set the longidutinal discretization */
  void setLongitudinalDisc( double zmin, double zmax, unsigned int Nz );

  /** Set SiO2 as the cladding material */
  void setCladdingSilicon( double energyInEv );

  /** Set a specific set of the cladding parametes. Not recommended. */
  void setCladdingDeltaBeta( double delta, double beta ); // This overrules the material properties

  /** Add alcohol inside */
  void setAlcoholInside( double energyInEv );

  /** Set ethylene glycol inside the waveguide */
  void setEthylenGlycolInside( double energyInEv );

  /** Add simualtion to save */
  void saveIndx( unsigned int indx );

  /** Set the signal length to be used in the FFT when computing the far field. Should be an 2^{some integer}*/
  void setFFTSignalLength( unsigned int length){fftSignalLength=length;};

  /** Add a visualizer to create a 2D colorplot of the wavefield */
  void setVisualizer( Visualizer &newvis ){ vis = &newvis; };

  /** Set directory to save the pictures. */
  void savePic( const char *dir );

  /** Do not generate random UID for the filename */
  void turnOffUID(){ generateUID = false; };

  /** Set the interval between each visualized wave field */
  void setDisplayInterval( unsigned int interval ){ displayEvery = interval; };
private:
  double getAngle( unsigned int indx ) const;
  PlaneWave pw;
  Cladding cladding;
  Cladding inside;
  ParaxialEquation eq;
  CrankNicholson solver;
  arma::mat farField;
  std::set<unsigned int> indxToSave;
  Visualizer *vis{NULL};
  std::string picDir;
  bool saveImages{false};
  bool vacuumInside{true};
  unsigned int displayEvery{10};
protected:
  /** Called during the solve process. Childs may want to process the far fields differently */
  virtual void processFarField(){};
  StraightWG2D wg;
  unsigned int fftSignalLength{65536};
  unsigned int nTheta{1};
  double thetaMin{0.0};
  double thetaMax{0.0};
  bool generateUID{true};
  unsigned int uid{0};

  /** Index in the FFT array corresponding to the exity angle in degree */
  unsigned int angleIndx( double angle ) const;

  /** Angle in degree corresponding to FFT index indx=(-N/2,N/2)*/
  double angleDeg( int indx ) const;
};
#endif
