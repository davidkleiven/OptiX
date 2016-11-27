#ifndef WAVE_GUIDE_FD_SIMULATION_H
#define WAVE_GUIDE_FD_SIMULATION_H
#include <cstddef>
#include <complex>
#include <string>
#include <jsoncpp/json/writer.h>
#include <armadillo>
#include "waveGuideBorder.hpp"
#include "borderTracker.hpp"
class Solver2D;

/** Struct storing the discretization parameters */
struct Disctretization
{
  double min;
  double max;
  double step;
};

typedef std::complex<double> cdouble;

class Cladding;
class ControlFile;
class WaveGuide1DSimulation;
class ParaxialSource;

/** Base class for 2D waveguide simulations */
class WaveGuideFDSimulation
{
public:
  /** Available boundaries in 2D */
  enum class Boundary_t {TOP, BOTTOM, LEFT, RIGHT};
  WaveGuideFDSimulation();
  WaveGuideFDSimulation(const char *name);
  ~WaveGuideFDSimulation();

  /** Set the transverse discretization */
  void setTransverseDiscretization( double xmin , double xmax, double step );

  /** Set the longitudinal discretization */
  void setLongitudinalDiscretization( double zmin , double zmax, double step );

  /** Get number of nodes in the transverse direction */
  unsigned int nodeNumberTransverse() const;

  /** Get number of nodes in the longitudinal direction */
  unsigned int nodeNumberLongitudinal() const;

  /** Get transverse discretiaztion */
  const Disctretization& transverseDiscretization() const{ return *xDisc; };

  /** Get longitudinal discretization */
  const Disctretization& longitudinalDiscretization() const { return *zDisc; };

  /** Compute far field */
  void computeFarField();

  /** Compute far field using a zero padded signal of a given length. Should be 2^{some integer} */
  void computeFarField( unsigned int signalLength ); // With padding to increase low freq resolution

  /** Get the wavenumber in nm^{-1}*/
  double getWavenumber() const{ return wavenumber; };

  /** Set wavenumber in nm^{-1}*/
  void setWavenumber( double k ){ wavenumber = k; };

  /** Set wavelength in nm */
  void setWaveLength( double lambda );

  /** Set 2D solver */
  void setSolver( Solver2D &solv );

  /** Set the length of the waveguide in nano meters */
  void setWaveguideLength( double wgl ){ wglength = wgl; };

  /** Get name of the waveguide simulation */
  std::string getName() const { return name; };

  /** Set cladding */
  void setCladding( const Cladding &clad );

  /** Set material inside the waveguide */
  void setInsideMaterial( const Cladding &mat );

  /** Get cladding */
  const Cladding& getCladding() const { return *cladding; };

  /** Run simulation */
  void solve();

  /** Save results to HDF5 files */
  void save( ControlFile &ctl ) const;

  /** Save results to HDF5 files, using only location where the intensity is above the threshold. Deprocated */
  void save( ControlFile &ctl, double intensityThreshold ) const;

  /** Save the waveguide borders */
  void saveWG( const std::string &fname ) const;

  /** Locate the waveguide borders */
  void extractWGBorders();

  /** Get intensity at position x,z */
  double getIntensity( double x, double z ) const; // Using linear interpolation

  /** Get intensity at array location ix, iz */
  double getIntensity( unsigned int ix, unsigned int iz ) const; // Returns value in matrix at (ix,iz)

  /** Get z-coordinate corresponding to the array index iz */
  double getZ( unsigned int iz ) const;

  /** Get x-coordinate corresponding to the array index ix */
  double getX ( int ix ) const;

  /** Get far field */
  const arma::vec& getFarField() const { return *farFieldModulus; };

  /** Get the exit field */
  void getExitField( arma::vec &vec ) const;

  /** Get 2D solver */
  const Solver2D& getSolver() const { return *solver; };

  /** Enable the use of border tracker */
  void useBorderTracker();

  /** Get a pointer to the border tracker */
  BorderTracker* getBorderTracker(){ return bTracker; };

  /** Get the array index closest to x, z */
  void closestIndex( double x, double z, unsigned int &ix, unsigned int &iz ) const;

  // Virtual methods
  /** Set incident field */
  virtual void setBoundaryConditions( const ParaxialSource& src ); // This function should fill the boundary

  /** Fill JSON object with parameters specific to this class */
  virtual void fillInfo( Json::Value &obj ) const {};

  /** Initialize. Relevant if loading an old solution */
  virtual void init( const ControlFile &ctl );

  /** Get the boundary condition at the specified boundary */
  virtual cdouble transverseBC( double z, Boundary_t bnd ) const;

  /** Get transverse boundary condition at position z. Can be used if is equal on x=xmin and x=xmax*/
  virtual cdouble transverseBC( double z ) const;

  /** Checks if the given point is inside the waveguide */
  virtual bool isInsideGuide( double x, double z ) const { return true; };
  // Refractive index: n = 1 - delta + i*beta
  /** Get the material properties */
  virtual void getXrayMatProp( double x, double z, double &delta, double &beta ) const;
protected:
  Solver2D *solver{NULL};
  Disctretization *xDisc; // Transverse
  Disctretization *zDisc; // Along optical axis
  double wavenumber;
  double wglength{1E10};
  std::string name;
  const Cladding *cladding{NULL};
  const Cladding *insideMaterial{NULL};
  arma::vec *farFieldModulus{NULL};
  BorderTracker *bTracker{NULL};

  /** Allocate solution matrix */
  double* allocateSolutionMatrix() const;

  /** Deallocate the solution matrix */
  void deallocateSolutionMatrix( double *matrix ) const;
  bool solverInitializedViaInit{false};

  /** Save only values whose intensity is above the threshold. Depricated */
  void sparseSave( const std::string &fname, double intensityThreshold ) const;

  /** Integrate from xStart to xEnd by a trapezoidal scheme */
  double trapezoidalIntegrateIntensityZ( unsigned int iz, unsigned int ixStart, unsigned int ixEnd ) const;

  /** Get exit field */
  void getExitField( arma::cx_vec &vec ) const;

  /** Save far field to HDF5 file */
  void saveFarField( const std::string &fname, unsigned int uid ) const;
  const ParaxialSource* src{NULL};
  std::vector<WaveGuideBorder> *wgborder{NULL};

  // Virtual funcitons
  /** Check if the point is after the waveguide end */
  virtual bool waveguideEnded( double x, double z ) const { return z > wglength; };
};
#endif
