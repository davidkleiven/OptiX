#ifndef POST_PROCESS_INTENSITY_H
#define POST_PROCESS_INTENSITY_H
#include "postProcessing.hpp"

class ParaxialSimulation;
namespace post
{

/** Module that gives the intensity of the solution */
class Intensity: public post::PostProcessingModule
{
public:
  Intensity(): post::PostProcessingModule("Intensity"){};

  /** Amplitude of the solution */
  void result( const Solver2D &solver, arma::mat &res ) override;
};

/** Module that returns the phase of the solution */
class Phase: public PostProcessingModule
{
public:
  Phase(): post::PostProcessingModule("Phase"){};

  /** Phase of the solution */
  void result( const Solver2D &solver, arma::mat &res ) override;
};

/** Module that computes the far field intensity pattern */
class FarField: public post::PostProcessingModule
{
public:
  FarField(): post::PostProcessingModule("FarField"){};

  /** Amplitude of the far field */
  void result( const Solver2D &solver, arma::vec &res ) override;

  /** Set pad length */
  void setPadLength( unsigned int pad ){ signalLength = pad; };

  /** Set the angular range */
  void setAngleRange( double angMin, double angMax );

  /** Link a paraxial simulation object */
  void linkParaxialSim( const ParaxialSimulation &simulation ){ sim = &simulation; };
private:
  unsigned int signalLength{0};
  double phiMin{-90.0};
  double phiMax{90.0};
  const ParaxialSimulation *sim{NULL};

  /** Computes the index in the far field array corresponding to a certain angle */
  unsigned int farFieldAngleToIndx( double angle, const arma::vec &res ) const;

  void reduceArray( arma::vec &res ) const;
};

};
#endif
