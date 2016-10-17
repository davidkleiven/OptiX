#ifndef WAVE_GUIDE_CURVED_2D_H
#define WAVE_GUIDE_CURVED_2D_H
#include "waveGuideFDSimulation.hpp"
#include <complex>
#include <vector>
#include <armadillo>

class ControlFile;

typedef std::complex<double> cdouble;

class CurvedWaveGuideFD: public WaveGuideFDSimulation
{
public:
  CurvedWaveGuideFD(): WaveGuideFDSimulation("CurvedWaveGuide2D"){};
  void getXrayMatProp( double x, double z, double &delta, double &beta ) const override final;
  void setBoundaryConditions() override final;
  void setRadiusOfCurvature( double newR ) { R = newR; };
  void setWidth( double newWidth ) { width = newWidth; };
  double getWidth() const { return width; };
  cdouble transverseBC( double z ) const override final;
  void computeTransmission( double step );
  void saveTransmission( ControlFile &ctl ) const;
  void getFieldInsideWG( arma::mat &matrix ) const;

  // Virtual functions
  virtual void fillInfo( Json::Value &obj ) const override;
  virtual void init( const ControlFile &ctl ) override;

  // Extracts the field inside the waveguide at distance wcrd from the edge
  // Thus: wcrd = 0.0: along the lower wall. wcrd=width extracts along the upper wall
  virtual void extractField( double wcrd, std::vector<cdouble> &res ) const {};
  //virtual void extractField( double wcrd, std::vector<double> &res ) const {} // real part of field
protected:
  CurvedWaveGuideFD( const char *name): WaveGuideFDSimulation(name){};
  double R;
  double width;
  virtual bool isInsideGuide( double x, double z ) const override;
  virtual double waveGuideStartX( double z ) const;
  virtual double waveGuideEndX( double z ) const;
  std::vector<double> transmission;
  double stepWhenComputingTransmission{0.0};
};
#endif
