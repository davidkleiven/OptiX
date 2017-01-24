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
  Intensity(): post::PostProcessingModule("amplitude", ReturnType_t::matrix2D){};

  /** Amplitude of the solution */
  virtual void result( const Solver &solver, arma::mat &res ) override;
};

/** Module that returns the phase of the solution */
class Phase: public PostProcessingModule
{
public:
  Phase(): post::PostProcessingModule("phase", ReturnType_t::matrix2D){};

  /** Phase of the solution */
  virtual void result( const Solver &solver, arma::mat &res ) override;
};

/** Module that computes the far field intensity pattern */
class FarField: public post::PostProcessingModule
{
public:
  FarField(): post::PostProcessingModule("farField", ReturnType_t::vector1D){};

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
class ExitField: public post::PostProcessingModule
{
public:
  ExitField():post::PostProcessingModule("exitField", ReturnType_t::vector1D){};

  /** Returns the real part of the exit field*/
  virtual void result( const Solver &solver, arma::vec &res ) override;
};

/** Module that returns the exit amplitude */
class ExitIntensity: public post::PostProcessingModule
{
public:
  ExitIntensity(): post::PostProcessingModule("exitIntensity", ReturnType_t::vector1D){};

  /** Returns the exit amplitude */
  virtual void result( const Solver &solver, arma::vec &res ) override;
};

/** Module for computing exit phase */
class ExitPhase: public post::PostProcessingModule
{
public:
  ExitPhase(): post::PostProcessingModule("exitPhase", ReturnType_t::vector1D){};

  /** Returns the phase at the exit */
  virtual void result( const Solver &solver, arma::vec &res ) override;
};

};
#endif
