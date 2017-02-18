#ifndef POST_PROCESS_INTENSITY_H
#define POST_PROCESS_INTENSITY_H
#include "postProcessing.hpp"
#include <armadillo>

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

  /** Amplitude of the solution 3D version */
  virtual void result( const Solver &solver, arma::cube &res ) override;
};

/** Module that returns the phase of the solution */
class Phase: public FieldQuantity
{
public:
  Phase(): post::FieldQuantity("phase"){};

  /** Phase of the solution */
  virtual void result( const Solver &solver, arma::mat &res ) override;

  /** Phase of the solution 3D version */
  virtual void result( const Solver &solver, arma::cube &res ) override;
};

/** Module that computes the far field intensity pattern */
class FarField: public post::ProjectionQuantity
{
public:
  enum class Pad_t{ZERO, LINEAR_FIT, CONSTANT};
  FarField(): post::ProjectionQuantity("farField"){};

  /** Set which type of padding to use */
  void setPadding( Pad_t padType ){ padding = padType; };

  /** Amplitude of the far field */
  virtual void result( const Solver &solver, arma::vec &res ) override;

  /** Far field in the 3D case */
  virtual void result( const Solver &solver, arma::mat &res ) override;

  /** Set pad length */
  void setPadLength( unsigned int pad ){ signalLength = pad; };

  /** Set reference solution to be subtracted off */
  void setReference( const arma::cx_mat &ref ){ reference = &ref; };

  /** Set the angular range */
  void setAngleRange( double angMin, double angMax );

  /** Link a paraxial simulation object */
  void linkParaxialSim( const ParaxialSimulation &simulation ){ sim = &simulation; };

  /** Add attributes */
  virtual void addAttrib( std::vector<H5Attr> &attr ) const override final;
private:
  enum class Dir_t{X,Y};
  unsigned int signalLength{0};
  double phiMin{-90.0};
  double phiMax{90.0};
  const ParaxialSimulation *sim{NULL};
  Pad_t padding{Pad_t::ZERO};
  const arma::cx_mat *reference{NULL};

  /** Computes the index in the far field array corresponding to a certain angle */
  unsigned int farFieldAngleToIndx( double angle, unsigned int size, Dir_t direction ) const;

  template<class T>
  void reduceArray( arma::Col<T> &res, Dir_t direction ) const;

  template<class T>
  void reduceArray( arma::Col<T> &res ) const;
  void reduceArray( arma::mat &res ) const;

  template<class T>
  void fftshift( arma::Col<T> &array );

  /** Pads the signal according to settings. Default is ZERO */
  void padSignal( arma::cx_vec &zeroPadded ) const;

  /** Linear padding */
  void linearPadding( arma::cx_vec &zeroPadded ) const;

  /** Pad with a constant on both sides */
  void constantPadding( arma::cx_vec &zeroPadded ) const;

  /** Subtract off the source */
  void subtractOffSource( arma::cx_mat &withoutSource ) const;
};

/** Module that returns the real part of the exit field */
class ExitField: public post::ProjectionQuantity
{
public:
  ExitField():post::ProjectionQuantity("exitField"){};

  /** Returns the real part of the exit field*/
  virtual void result( const Solver &solver, arma::vec &res ) override;

  /** Exit field, 3D version */
  virtual void result( const Solver &solver, arma::mat &res ) override;
};

/** Module that returns the exit amplitude */
class ExitIntensity: public post::ProjectionQuantity
{
public:
  ExitIntensity(): post::ProjectionQuantity("exitIntensity"){};

  /** Returns the exit amplitude */
  virtual void result( const Solver &solver, arma::vec &res ) override;

  /** Exit amplitude 3D version */
  virtual void result( const Solver &solver, arma::mat &res ) override;
};

/** Module for computing exit phase */
class ExitPhase: public post::ProjectionQuantity
{
public:
  ExitPhase(): post::ProjectionQuantity("exitPhase"){};

  /** Returns the phase at the exit */
  virtual void result( const Solver &solver, arma::vec &res ) override;

  /** Returns the phase at the exit, 3D version */
  virtual void result( const Solver &solver, arma::mat &res ) override;
};

};
#endif
