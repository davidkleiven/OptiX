#ifndef SOLVER2D_H
#define SOLVER2D_H
#include <string>
#include <complex>
#include <json/writer.h>
#include <armadillo>
#include <visa/gaussianKernel.hpp>
#include <visa/lowPassFilter.hpp>

class WaveGuideFDSimulation;
class ParaxialEquation;
class ParaxialSimulation;
class BoundaryCondition;

typedef std::complex<double> cdouble;

/** Base class for the 2D solvers */
class Solver2D
{
public:
  /** Enum for real or imaginary components */
  enum class Comp_t{REAL,IMAG};
  Solver2D( const char* name ):name(name){};
  virtual ~Solver2D();

  /** Get the name of the solution */
  std::string getName() const { return name; };

  /** Set waveguide to operate on */
  void setSimulator( ParaxialSimulation &guide );

  const ParaxialSimulation& getSimulator() const { return *guide; };

  /** Set paraxial equation to solve */
  void setEquation( const ParaxialEquation &equation ){ eq = &equation; };

  /** Get the solution. Depricated */
  const arma::cx_mat& getSolution( unsigned int iz ) const { return *solution; }; // Depricated. iz is not used.

  /** Get the solution */
  const arma::cx_mat& getSolution() const { return *solution; };

  /** Import solution from HDF5 file */
  bool importHDF5( const std::string &fname );

  /** Import solution from HDF5 files containing the real part and the amplitude. Depricated */
  bool importHDF5( const std::string &realpart, const std::string &amplitude );

  /** Get the real part of the solution */
  void realPart( double *realsolution ) const;

  /** Get the imaginary part of the solution */
  void imagPart( double *imagsolution ) const;

  /** Get the solution in matrix form */
  void getField( arma::mat &field ) const;

  /** Get the phase of the solution in matrix */
  void getPhase( arma::mat &phase ) const;

  // Set boundary condition at z=0
  /** Set boundary condition at z = 0*/
  void setLeftBC( const cdouble values[] );

  /** Set boundary condition at x=xmin and x=xmax */
  void setXBC( const cdouble valuesTop[], const cdouble valuesBottom[] ); // BC at top (x=xmax) and bottom (x=xmin)

  /** Get last solution. NOTE: Return prev solution as this is not filtered */
  arma::cx_vec& getLastSolution() const { return *prevSolution; };

  /** Run simulation */
  void solve();

  /** Perform one step */
  void step();

  /** Filter and downsample in the longitudinal direction */
  void filterInLongitudinalDirection();

  /** Down sample in longitudinal direction */
  void downSampleLongitudinalDirection();

  /** Reset the counter */
  void reset(){ currentStep = 1; };

  /** Sets the transverse boundary conditions */
  void addBoundaryCondition( const BoundaryCondition &bCond ){ bc = &bCond; };

  // Virtual functions
  /** Pure virtual function for solving the system */

  /** Fill JSON object with parameters specific to this class */
  virtual void fillInfo( Json::Value &obj ) const;
protected:
  std::string name;
  ParaxialSimulation *guide;
  const ParaxialEquation *eq{NULL};
  arma::cx_mat *solution{NULL};
  arma::cx_vec *prevSolution{NULL};
  arma::cx_vec *currentSolution{NULL};
  const BoundaryCondition *bc{NULL};
  unsigned int currentStep{1};

  unsigned int Nx{0};
  unsigned int Nz{0};
  double stepX{1.0}, stepZ{1.0};
  double xmin{0.0};
  double zmin{0.0};
  double wavenumber{1.0};

  /** Get the solution */
  arma::cx_mat& getSolution( unsigned int iz ) { return *solution; };

  /** Get solution, rear or imaginary part specified by comp */
  void realOrImagPart( double *solution, Comp_t comp ) const;

  /** Copies the solution in to the solution matrix */
  void copyCurrentSolution( unsigned int step );

  /** Solve one step */
  virtual void solveStep( unsigned int step ) = 0;

  /** Set the required parameters from the waveguide object */
  void initValuesFromWaveGuide();

  visa::GaussianKernel kernel;
  visa::LowPassFilter filter;
};
#endif
