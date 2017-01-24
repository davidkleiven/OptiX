#ifndef POST_PROCESS_INTENSITY_H
#define POST_PROCESS_INTENSITY_H
#include "postProcessing.hpp"

class ParaxialSimulation;
namespace post
{

/** Module that gives the intensity of the solution */
class Intensity: public post::FieldQuantity
{
public:
  Intensity(): post::FieldQuantity("amplitude"){};

  /** Amplitude of the solution */
  virtual void result( const Solver &solver, arma::mat &res ) override;
};

/** Module that returns the phase of the solution */
class Phase: public FieldQuantity
{
public:
  Phase(): post::FieldQuantity("phase"){};

  /** Phase of the solution */
  virtual void result( const Solver &solver, arma::mat &res ) override;
};

/** Module that computes the far field intensity pattern */
class FarField: public post::ProjectionQuantity
{
public:
  FarField(): post::ProjectionQuantity("farField"){};

  /** Amplitude of the far field */
  virtual void result( const Solver &solver, arma::vec &res ) override;

  /** Set pad length */
  void setPadLength( unsigned int pad ){ signalLength = pad; };

  /** Set the angular range */
  void setAngleRange( double angMin, double angMax );

  /** Link a paraxial simulation object */
  void linkParaxialSim( const ParaxialSimulation &simulation ){ sim = &simulation; };

  /** Add attributes */
  virtual void addAttrib( std::vector<H5Attr> &attr ) const override final;
private:
  unsigned int signalLength{0};
  double phiMin{-90.0};
  double phiMax{90.0};
  const ParaxialSimulation *sim{NULL};

  /** Computes the index in the far field array corresponding to a certain angle */
  unsigned int farFieldAngleToIndx( double angle, const arma::vec &res ) const;

  void reduceArray( arma::vec &res ) const;
};

/** Module that returns the real part of the exit field */
class ExitField: public post::ProjectionQuantity
{
public:
  ExitField():post::ProjectionQuantity("exitField"){};

  /** Returns the real part of the exit field*/
  virtual void result( const Solver &solver, arma::vec &res ) override;
};

/** Module that returns the exit amplitude */
class ExitIntensity: public post::ProjectionQuantity
{
public:
  ExitIntensity(): post::ProjectionQuantity("exitIntensity"){};

  /** Returns the exit amplitude */
  virtual void result( const Solver &solver, arma::vec &res ) override;
};

/** Module for computing exit phase */
class ExitPhase: public post::ProjectionQuantity
{
public:
  ExitPhase(): post::ProjectionQuantity("exitPhase"){};

  /** Returns the phase at the exit */
  virtual void result( const Solver &solver, arma::vec &res ) override;
};

};
#endif
