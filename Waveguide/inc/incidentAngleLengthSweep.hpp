#ifndef INCIDENT_ANGLE_LENGTH_SWEEP_H
#define INCIDENT_ANGLE_LENGTH_SWEEP_H
#include "incidentAngleSweep.hpp"
#include <vector>
#include <armadillo>
#include <string>

/** Stores the far field obtained at different positions */
struct FarFieldDataset
{
  double length;
  arma::mat farField;
};

/** Class that computes the far field at several waveguide lengths */
class IncidentAngleLengthSweep: public IncidentAngleSweep
{
public:
  IncidentAngleLengthSweep(){};

  /** Set minium length to compute far field (maximum is given by incident angle's waveguide length)*/
  void setLmin( double lmin ){ Lmin=lmin; nextZ=lmin; };

  /** Set number of lengths for which the far field should be computed */
  void setNumberOfLengths( unsigned int N ){ Nlengths=N; };

  /** Save all the datasets to a file */
  virtual void save( const std::string &fname ) const override;
private:
  double Lmin{2.0};
  double Nlengths{100};
  unsigned int currentLengthIter{0};
  std::vector<FarFieldDataset> data;
  unsigned int currentTheta{0};
  double nextZ;

  /** Compute far field at base on the field from several locations */
  virtual void processStep( double z ) override;

  /** Called when the solver proceeds to the next angle */
  virtual void proceedToNextAngle() override;
};
#endif
