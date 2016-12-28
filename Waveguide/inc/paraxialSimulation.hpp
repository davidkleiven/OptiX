#ifndef PARAXIAL_SIMULATION_H
#define PARAXIAL_SIMULATION_H
#include <complex>
#include <H5Cpp.h>
#include <armadillo>
#include <json/writer.h>
#include "farFieldParameters.hpp"
#include "h5Attribute.hpp"
#include <vector>
class Solver2D;
class ControlFile;
class ParaxialSource;
class BorderTracker;

/** Struct storing the discretization parameters */
struct Disctretization
{
  double min;
  double max;
  double step;
  unsigned int downsamplingRatio{1};
};
typedef std::complex<double> cdouble;
/** Base class for all paraxial simulations */
class ParaxialSimulation
{
public:
  /** Available boundaries in 2D */
  enum class Boundary_t {TOP, BOTTOM, LEFT, RIGHT};
  explicit ParaxialSimulation( const char* name );
  virtual ~ParaxialSimulation();

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

  /** Set the padding value used when computing the far field */
  void setFarFieldPadValue( double padValue ){ farParam.padValue=padValue; };

  /** Set the angle range where the far field will be stored */
  void setFarFieldAngleRange( double phiMin, double phiMax );

  /** Compute far field */
  void computeFarField();

  /** Compute far field using a zero padded signal of a given length. Should be 2^{some integer} */
  void computeFarField( unsigned int signalLength ); // With padding to increase low freq resolution

  /** Compute far field based on the field at z-position given by pos */
  void computeFarField( unsigned int signalLength, double pos );

  /** Get the wavenumber in nm^{-1}*/
  double getWavenumber() const{ return wavenumber; };

  /** Get the wavelength in nm */
  double getWavelength() const;

  /** Get wave energy in eV */
  double getEnergy() const;

  /** Set wavenumber in nm^{-1}*/
  void setWavenumber( double k ){ wavenumber = k; };

  /** Set wavelength in nm */
  void setWaveLength( double lambda );

  /** Set 2D solver */
  void setSolver( Solver2D &solv );

  /** Get name of the waveguide simulation */
  std::string getName() const { return name; };

  /** Run simulation */
  void solve();

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

  /** Get the array index closest to x, z */
  void closestIndex( double x, double z, unsigned int &ix, unsigned int &iz ) const;

  /** Enable/disable storing of the intensity and phase for a contour plot */
  void saveContour( bool save=true ){ saveColorPlot=save; };

  // Virtual methods
  /** Set incident field */
  virtual void setBoundaryConditions( const ParaxialSource& src ); // This function should fill the boundary

  /** Fill JSON object with parameters specific to this class */
  virtual void fillInfo( Json::Value &obj ) const {};

  /** Initialize. Relevant if loading an old solution */
  virtual void init( const ControlFile &ctl ){}; // TODO: Implement this

  /** Get the boundary condition at the specified boundary */
  virtual cdouble transverseBC( double z, Boundary_t bnd ) const{};

  /** Get transverse boundary condition at position z. Can be used if is equal on x=xmin and x=xmax*/
  virtual cdouble transverseBC( double z ) const{};

  /** Get the material properties */
  virtual void getXrayMatProp( double x, double z, double &delta, double &beta ) const{};

  /** Save results to HDF5 files */
  virtual void save( ControlFile &ctl );

  /** Return a border tracker object. Only relevant for geometries that tracks the border i.e. waveguides */
  virtual BorderTracker* getBorderTracker(){ return NULL; };
protected:
  Solver2D *solver{NULL};
  Disctretization *xDisc; // Transverse
  Disctretization *zDisc; // Along optical axis
  arma::vec *farFieldModulus{NULL};
  double wavenumber;
  std::string name;
  const ParaxialSource* src{NULL};
  H5::H5File* file{NULL};
  std::vector<std::string> dsetnames;
  bool solverInitializedViaInit{false};
  bool saveColorPlot{true};
  FarFieldParameters farParam;
  std::vector<H5Attr> commonAttributes;

  /** Get exit field */
  void getExitField( arma::cx_vec &vec ) const;

  /** Save far field to HDF5 file */
  void saveFarField();

  /** Extra calling */
  virtual void saveSpecialDatasets( hid_t file_id, std::vector<std::string> &dset ) const{};

  /** Add armadillo matrix to HDF5 file */
  void saveArmaMat( const arma::mat &matrix, const char* dsetname, const std::vector<H5Attr> &attr );

  /** Add armadillo matrix to HDF5 using the common attributes */
  void saveArmaMat( const arma::mat &matrix, const char* dsetname );

  /** Add armadillo vector to HDF5 */
  void saveArmaVec( const arma::vec &vec, const char* dsetname, const std::vector<H5Attr> &attr );

  /** Add armadillo vector to HDF5 using the default attributes */
  void saveArmaVec( const arma::vec &vec, const char* dsetname );

  /** Add STL vector to HDF5 */
  void saveVec( const std::vector<double> &vec, const char* dsetname );

  /** Add attribute to dataset */
  void addAttribute( H5::DataSet &ds, const char* name, double value );

  /** Add attribute to dataset */
  void addAttribute( H5::DataSet &ds, const char* name, int value );

  /** Extracts the part of the far field corresponding to the angles in far field parameters */
  void extractFarField( arma::vec &newFarField ) const;

  /** Computes the index in the far field array corresponding to a certain angle */
  unsigned int farFieldAngleToIndx( double angle ) const;

  /** Pad the exit signal */
  virtual cdouble padExitField( double x, double z ) const { return farParam.padValue; };
};
#endif
