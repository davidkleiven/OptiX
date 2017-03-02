#ifndef GENERIC_SCATTERING_H
#define GENERIC_SCATTERING_H
#include "paraxialSimulation.hpp"
#include "postProcessMod.hpp"
#include "fftSolver3D.hpp"
#include "gaussianBeam.hpp"
#include "materialFunction.hpp"

class GenericScattering: public ParaxialSimulation
{
public:
  GenericScattering( const char* name );

  void setMaterial(  const MaterialFunction &mat ){ material = &mat; };

  /** Set the waist of the Gaussian beam */
  void setBeamWaist( double waist ){ gbeam.setWaist(waist); };

  /** Sets the maximum scattering angle that will be stored */
  void setMaxScatteringAngle( double angMax );

  /** Solve scattering */
  virtual void solve() override;

  /** Gets the X-ray material properties */
  virtual void getXrayMatProp( double x, double y, double z, double &delta, double &beta ) const override;

  /** Prints all the attributes */
  void printInfo() const;

  // Attributes
  unsigned int exportNx{512};
  unsigned int exportNy{512};
  unsigned int FFTPadLength{1024};
  double xmin{-1.0}, xmax{1.0}, ymin{-1.0}, ymax{1.0}, zmin{-1.0}, zmax{1.0};
  double dx{0.1}, dy{0.1}, dz{0.1};
  unsigned int downSampleX{1}, downSampleY{1}, downSampleZ{1};
  double wavelength{0.1569};
  double intensityMin{0.0}, intensityMax{1.0};
  double phaseMin{-3.14159}, phaseMax{3.14159};
  std::string imgname{""};
  bool realTimeVisualization{false};
private:
  const MaterialFunction *material{NULL};
  post::ExitField ef;
  post::ExitIntensity ei;
  post::ExitPhase ep;
  post::FarField ff;

  GaussianBeam gbeam;
  FFTSolver3D fft3Dsolver;

  arma::cx_mat reference;
  bool isReferenceRun{true};

  /** Has to be called before simulation is solved */
  void init();
};
#endif
