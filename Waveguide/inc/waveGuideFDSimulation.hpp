#ifndef WAVE_GUIDE_FD_SIMULATION_H
#define WAVE_GUIDE_FD_SIMULATION_H
#include <cstddef>
#include <complex>
#include <string>
#include <json/writer.h>
#include <armadillo>
#include "waveGuideBorder.hpp"
#include "borderTracker.hpp"
#include "paraxialSimulation.hpp"
class Solver2D;
typedef std::complex<double> cdouble;

class Cladding;
class ControlFile;
class WaveGuide1DSimulation;
class ParaxialSource;

/** Base class for 2D waveguide simulations */
class WaveGuideFDSimulation: public ParaxialSimulation
{
public:
  WaveGuideFDSimulation();
  WaveGuideFDSimulation(const char *name);
  virtual ~WaveGuideFDSimulation();

  /** Set the length of the waveguide in nano meters */
  void setWaveguideLength( double wgl ){ wglength = wgl; };

  /** Set cladding */
  void setCladding( const Cladding &clad );

  /** Set material inside the waveguide */
  void setInsideMaterial( const Cladding &mat );

  /** Get cladding */
  const Cladding& getCladding() const { return *cladding; };

  /** Save the waveguide borders */
  void saveWG( const std::string &fname ) const;

  /** Locate the waveguide borders */
  void extractWGBorders();

  /** Enable the use of border tracker */
  void useBorderTracker();

  /** Get a pointer to the border tracker */
  BorderTracker* getBorderTracker() override { return bTracker; };

  /** Get pointer to const border tracker */
  const BorderTracker* getBorderTracker() const { return bTracker; };

  // Virtual methods

  /** Fill JSON object with parameters specific to this class */
  virtual void fillInfo( Json::Value &obj ) const override{};

  /** Initialize. Relevant if loading an old solution */
  virtual void init( const ControlFile &ctl ) override{}; // TODO: Implement this

  /** Get the boundary condition at the specified boundary */
  virtual cdouble transverseBC( double z, Boundary_t bnd ) const override;

  /** Get transverse boundary condition at position z. Can be used if is equal on x=xmin and x=xmax*/
  virtual cdouble transverseBC( double z ) const override;

  /** Checks if the given point is inside the waveguide */
  virtual bool isInsideGuide( double x, double z ) const { return true; };
  // Refractive index: n = 1 - delta + i*beta
  /** Get the material properties */
  virtual void getXrayMatProp( double x, double z, double &delta, double &beta ) const override;

  /** Save datasets specific to this class */
  virtual void saveSpecialDatasets( hid_t file_id, std::vector<std::string> &dset ) const override;

  /** Save datasets */
  virtual void save( ControlFile &ctl ) override;
protected:
  double wglength{1E10};
  const Cladding *cladding{NULL};
  const Cladding *insideMaterial{NULL};
  BorderTracker *bTracker{NULL};

  /** Allocate solution matrix */
  double* allocateSolutionMatrix() const;

  /** Deallocate the solution matrix */
  void deallocateSolutionMatrix( double *matrix ) const;

  /** Save only values whose intensity is above the threshold. Depricated */
  void sparseSave( const std::string &fname, double intensityThreshold ) const;

  /** Integrate from xStart to xEnd by a trapezoidal scheme */
  double trapezoidalIntegrateIntensityZ( unsigned int iz, unsigned int ixStart, unsigned int ixEnd ) const;

  std::vector<WaveGuideBorder> *wgborder{NULL};

  // Virtual funcitons
  /** Check if the point is after the waveguide end */
  virtual bool waveguideEnded( double x, double z ) const { return z > wglength; };

  /** Pad the signal corresponding to the unaffected source */
  virtual cdouble padExitField( double x, double z ) const override;
};
#endif
