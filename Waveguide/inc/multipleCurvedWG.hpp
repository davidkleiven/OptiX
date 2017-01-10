#ifndef MULTIPLE_CURVED_WG_H
#define MULTIPLE_CURVED_WG_H
#include "paraxialSimulation.hpp"
#include "paraxialEquation.hpp"
#include "paraxialEqNegative.hpp"
#include "crankNicholson.hpp"
#include "planeWave.hpp"
#include "cladding.hpp"
#include "fixedValuesSource.hpp"
#include <string>
#include <vector>
#include <map>

class CurvedWGConfMap;
class MultipleCurvedWG: public ParaxialSimulation
{
public:
  MultipleCurvedWG(): ParaxialSimulation("MultipleCurvedWG"){};
  ~MultipleCurvedWG();

  /** Loads the geometry from a JSON file */
  void loadWaveguides( const std::string &jsonfname );

  /** Initializes the simulation */
  void init( const std::map<std::string,double> &params );

  /** Returns the intensity */
  const arma::mat& getIntensity() const { return *intensity; };

  /** Solve the system */
  virtual void solve() override;
private:
  std::vector<CurvedWGConfMap*> *waveguides{NULL};
  std::vector<double> angles;
  ParaxialEquation eq;
  ParaxialEqNegative eqNeg;
  CrankNicholson solver;
  PlaneWave pw;
  FixedValuesSource fsource;
  Cladding cladding;
  arma::mat *intensity;
  arma::vec *transmittivity;
  unsigned int NzNextFillStartIntensity{0};
  unsigned int NzNextFillStartTrans{0};

  template<class elemType>
  void flipWrtCenterOfWG( elemType vec[], unsigned int N ) const;

  void flipWrtCenterOfWG( arma::cx_vec &vec ) const;

  void flipWrtCenterOfWG( arma::mat &mat ) const;
  //void flipWrtCenterOfWG( arma::mat &mat ) const;

  void processSolution( CurvedWGConfMap &wg );
};
#endif
