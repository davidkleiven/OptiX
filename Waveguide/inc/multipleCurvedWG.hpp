#ifndef MULTIPLE_CURVED_WG_H
#define MULTIPLE_CURVED_WG_H
#include "paraxialSimulation.hpp"
#include "paraxialEquation.hpp"
#include "crankNicholson.hpp"
#include "planeWave.hpp"
#include "cladding.hpp"
#include "arraySource.hpp"
#include "controlFile.hpp"
#include "postProcessMod.hpp"
#include <string>
#include <vector>
#include <map>
#include <json/reader.h>

class CurvedWGConfMap;
class MultipleCurvedWG: public ParaxialSimulation
{
public:
  MultipleCurvedWG(): ParaxialSimulation("MultipleCurvedWG"){};
  ~MultipleCurvedWG();

  /** Loads the geometry from a JSON file */
  void loadWaveguides( const std::string &jsonfname );

  /** Initialize the waveguides from a JSON object */
  void loadWaveguides( const Json::Value &root );

  /** Initializes the simulation */
  void init( const std::map<std::string,double> &params );

  /** Returns the intensity */
  const arma::mat& getIntensity() const { return *intensity; };

  /** Get the transmittivity */
  const arma::vec& getTransmittivity() const { return *transmittivity; };

  /** Solve the system */
  virtual void solve() override;

  /** Save run */
  virtual void save( ControlFile &ctl ) override;
private:
  std::vector<CurvedWGConfMap*> *waveguides{NULL};
  std::vector<double> angles;
  ParaxialEquation eq;
  ArraySource fsource;
  Cladding cladding;
  arma::mat *intensity;
  arma::vec *transmittivity;
  unsigned int NzNextFillStartIntensity{0};
  unsigned int NzNextFillStartTrans{0};
  std::string imagefile;
  std::string geometryfile;
  post::FarField farfield;
  post::ExitField exitfield;
  post::ExitPhase exPhase;

  /** Flips array with respect to the center of the waveguide */
  template<class elemType>
  void flipWrtCenterOfWG( elemType vec[], unsigned int N ) const;

  /** Flips Armadillo vector with respect to the center of the waveguide */
  void flipWrtCenterOfWG( arma::cx_vec &vec ) const;

  /** Flips Armadillo matrix with respect to the center of the waveguide */
  void flipWrtCenterOfWG( arma::mat &mat ) const;

  /** Computes the required phase difference between two waveguides */
  double phaseDifference( const CurvedWGConfMap &source, const CurvedWGConfMap &target ) const;

  void processSolution( CurvedWGConfMap &wg );
};
#endif
