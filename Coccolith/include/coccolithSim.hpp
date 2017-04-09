#ifndef COCCOLITH_SIMULATION_H
#define COCCOLITH_SIMULATION_H
#include "config.h"
#include <string>
#include "voxelMaterial.hpp"
#include "fieldMonitors.hpp"
#include <complex>
#include <json/writer.h>
#include "stokesParameters.hpp"
#include <array>
#ifdef HAVE_LIB_VISA
  #include <visa/visa.hpp>
#endif
#include <cstdlib>
#define UID_MAX 10000000

typedef std::complex<double> cdouble;
enum class MainPropDirection_t{ X, Y, Z };
typedef MainPropDirection_t IntegrationDir_t;
typedef MainPropDirection_t RotationAxis_t;
enum class Plane_t{XY, XZ, YZ};

enum class SourcePosition_t{TOP, BOTTOM};

/** Struct for holding the local Stokes parameters */
struct LocalStokes
{
  std::vector<double> I;
  std::vector<double> Q;
  std::vector<double> U;
  std::vector<double> V;

  // Store a value for each of the electric field components as well
  double Ephi;
  double Etheta;
};

class Stokes
{
public:
  Stokes(){};
  Stokes(int I, int Q, int U, int V):I(I),Q(Q),U(U),V(V){};
  Stokes( int vec[4] ):I(vec[0]), Q(vec[1]), U(vec[2]), V(vec[3]){};

  /** Rotate the stokes vector by angle in radians */
  void rotate( double angleRad );

  bool operator==(const Stokes &other) const;
  double I{1};
  double Q{1};
  double U{0};
  double V{0};
};

class CoccolithSimulation
{
public:
  CoccolithSimulation(){};
  virtual ~CoccolithSimulation();

  /** LSet voxel material to use in the simulation */
  void setMaterial( VoxelMaterial &mat ){ material = &mat; };

  /** Set the incident wavevector. Note that this is a static member */
  void setIncWaveVector( const meep::vec &wave ){ waveVec = wave; };

  /** Sets the initial stokes vector */
  void setIncStokesVector( const int stk[4] );

  /** Set the main propagation direction */
  void setMainPropagationDirection( MainPropDirection_t propDir );

  /** Return the wavelength in pixel units */
  double getWavelength() const;

  /** Set the PML thickness in number of wavelengths */
  void setPMLInWavelengths( double newThick );

  /** Run the simulation */
  virtual void run();

  /** Initialize the source profile */
  void initSource( double freq, double fwidth );

  /** Initialize the simulation */
  virtual void init();

  /** Set the number of frequencies to be used in the harmonic inversion */
  void setNfreqFT( unsigned int numberOfFreq ){ nfreq = numberOfFreq; };

  /** Set the number of timesteps to save */
  void setNumberOfImgSave( unsigned int newnSave){ nSave = newnSave; };

  /** Set side for the source. TOP: Close to the border indexed 0, BOTTOM: close to the opposite border */
  void setSourceSide( SourcePosition_t newPos ){ srcPos = newPos; };

  /** Set the end time */
  void setEndTime( double t );

  /** Print info about the domain */
  void domainInfo() const;

  /** Exporting results */
  virtual void exportResults();

  /** The run is a reference run */
  void runWithoutScatterer();

  /** Add scatterer */
  void runWithScatterer();

  /** Resets the simulation */
  void reset();

  /** Sets the frequency between every time the graphics is updated */
  void setPlotUpdateFreq( unsigned int everyIter ){ plotUpdateFreq = everyIter; };

  /** Runs without visualization */
  void disableRealTimeVisualization(){ realTimeVisualization = false; };

  /** Computes the scattering assymetry factor for each frequency based on fields at distance R */
  void scatteringAssymmetryFactor( std::vector<double> &g, double R, unsigned int Nsteps );

  /** Computes and store the far field quantities */
  void farFieldQuantities();

  /** Loads previously stored near field currents */
  void loadBoundingCurrents( const char* fname );

  /** Evaluate far field */
  //void farFieldOnBox;

  /** Returns a reference to the structure */
  meep::structure& getStructure(){ return *struc; };

  double resolution{1.0};
  std::string uid{""};
  double additionalVaccumLayerPx{0.0}; // Additional vacuum layer that will be added outside the region

  /** Set a Sellmeier material */
  void setSellmeierMaterial( const SellmeierMaterial &mat ){ sellmeier = &mat; };

  /** Prefix that will be added to the filename */
  std::string prefix{""};

  /** Set to true if the scattering assummetry factor should be computed */
  bool computeAsymmetryFactor{false};

  /** Use an 45 degree polarization in the sources */
  bool computeStokesParameters{false};

  /** Set the order of the Gauss Legendre that is used for the asymmetry function calculation */
  unsigned int gaussLegendreOrder{17};

  /** Number of azimuthal angles to average over */
  unsigned int numberOfAzimuthalSteps{3};
protected:
  static Stokes incStoke;
  bool initialStokesVectorSet{false};
  static std::array<Stokes,6> supportedStokes;
  VoxelMaterial *material{NULL};
  const SellmeierMaterial *sellmeier{NULL};
  MainPropDirection_t propagationDir{MainPropDirection_t::Z};
  meep::component fieldComp{meep::Ex};
  meep::component secondComp{meep::Ey};

  meep::volume* srcVol{NULL};
  meep::structure* struc{NULL};
  meep::grid_volume gdvol;
  meep::fields* field{NULL};
  meep::src_time *sourceTime{NULL}; // Deleted by MEEP
  meep::gaussian_src_time *source{NULL};
  std::string outdir{"data"};
  unsigned int nSave{30};
  bool isInitialized{false};
  unsigned int plotUpdateFreq{30};
  Json::Value root;

  double pmlThicknessInWavelengths{3.0};
  double centerFrequency{1.0};
  double freqWidth{0.5};
  bool materialLoaded{false};
  unsigned int nfreq{100};
  unsigned int nMonitorX{256};
  unsigned int nMonitorY{256};
  unsigned int nMonitorZ{256};
  double tEnd{100.0};
  #ifdef HAVE_LIB_VISA
    visa::WindowHandler *plots{NULL};
  #endif
  bool userOverridedEndTime{false};
  bool realTimeVisualization{true};
  bool geoIsInitialized{false};
  std::string reflFluxPlaneBackup{"reflectedFlux"};
  std::string reflFluxBoxBackup{"reflectedFluxBox"};
  std::string n2fBoxBackup{"n2fBox"};

  /** Visualized intensity */
  void visualize();

  // Flux planes
  meep::volume *dftVolTransmit{NULL};
  meep::volume *dftVolRefl{NULL};
  meep::volume *dftVolBox{NULL};
  meep::volume_list *faces{NULL}; // N2F faces
  meep::dft_flux *transmitFlux{NULL};
  meep::dft_flux *reflFlux{NULL};
  meep::dft_flux *fluxBox{NULL};
  meep::dft_near2far *n2fBox{NULL};

  // Monitor planes 1
  FieldMonitor *monitor1{NULL};
  FieldMonitor *monitor2{NULL};

  arma::mat bkg1;
  arma::mat bkg2;
  arma::mat radialPoyntingVector;
  std::vector<arma::mat> stokesI;
  std::vector<arma::mat> stokesQ;
  std::vector<arma::mat> stokesV;
  std::vector<arma::mat> stokesU;
  std::vector<double> stokesIInc;
  std::vector<double> stokesQInc;
  std::vector<double> stokesUInc;
  std::vector<double> stokesVInc;
  arma::mat Ephi;
  arma::mat Etheta;
  arma::mat stokesIAzim;
  arma::mat stokesQAzim;
  arma::mat stokesUAzim;
  arma::mat stokesVAzim;
  unsigned int currentTheta{0};
  std::vector<double> thetaValues;
  std::vector<double> phiValues;

  // Source corners
  meep::vec crn1;
  meep::vec crn2;

  SourcePosition_t srcPos{SourcePosition_t::TOP};

  meep::h5file *file{NULL};

  unsigned int totalNumberOfFarFieldEvaluations{0};
  unsigned int numberOfTimesThereIsRadialFieldComponent{0};
  double relativeRadFieldCompThreshold{1E-2};
  double currentStokesVectorRotationAngleRad{0.0};

  /** Add a source volume */
  void addSourceVolume();

  /** Add a source to the field */
  void addSource();

  /** Add a structure for the simulation */
  void addStructure();

  /** Add fields */
  void addFields();

  /** Add DFT flux planes */
  void addFluxPlanes();

  /** Get the position of the source plane */
  double getSrcPos() const;

  /** Get the position of the flux plane next to the source */
  double getSrcFluxPos() const;

  /** Get the position of the flux plane after the scatterer */
  double getTransFluxPos() const;

  /** Returns the position of the lower indexed border in the propagation direction */
  double getLowerBorderInPropDir() const;

  /** Returns the position of the upper indexed border in the propagation direction */
  double getUpperBorderInPropDir() const;

  /** Returns the PML thickness in pixel units */
  double getPMLThickness() const;

  /** Sets two monitor planes passing through the center of the computational domain */
  void setMonitorPlanes();

  /** Project the epsilon along specified axis */
  void projectedEpsilon( arma::mat &matrix, IntegrationDir_t dir );

  /** Stores the results from the DFT spectrums */
  void saveDFTSpectrum();

  /** Save DFT parameters */
  void saveDFTParameters();

  /** Save parameters specific to the geometry and source */
  void saveGeometry();

  /** Computes the position corresponding to position in matrix given a certin plane*/
  void getPos( unsigned int row, unsigned int col, Plane_t proj, double dx, double dy, double dz, meep::vec &pos ) const;

  /** Returns an estimated time to cross the domain in MEEP units */
  double estimatedTimeToPropagateAcrossDomain() const;

  /** Initialize the geometry based on */
  void initializeGeometry();

  static meep::vec waveVec;

  static cdouble amplitude( const meep::vec &r );
  static cdouble amplitude2( const meep::vec &r );

  /** Sets the UID based on local time */
  void setUID();

  /** Adds susceptibilities to the structure */
  void updateStructure();

  /** Adds near to far field planes */
  void addN2FPlanes( const meep::volume &box );

  /** Integrate the radial component of the poyntings vector in the azimuthal direction */
  void azimuthalIntagration( double R, double theta, unsigned int Nsteps, std::vector<double> &res );

  /** Permute x,y,z given with propagation direction in Z to fit other propgation directions */
  template<class T>
  void permumteToFitPropDir( T &x, T &y, T &z ) const;

  /** Updates the Stokes parameter array based on the fields, weight is the Legendre-Gauss integration weight */
  void updateStokesParameters( const cdouble EH[], unsigned int evalPointIndx, double weight );

  /** Computes the contribution to the phase function */
  double phaseFunctionContribution( const cdouble[3] ) const;

  /** Redfine theta be the compliment based on the propagation direction */
  bool redefineTheta() const;

  /** Computes the two E-hat vectors orthonormal to the propagation direction given by theta, phi (in radians) */
  void computeEvectorOrthogonalToPropagation( double theta, double phi, meep::vec &E1hat, meep::vec &E2hat ) const;

  /** Computes the two E-hat vectors orthonormal to the propagation direction */
  void computeEvectorOrthogonalToPropagation( const meep::vec &r, meep::vec &E1hat, meep::vec &E2hat );

  /** Computes the Stokes parameters in the direction given by theta and phi */
  void getLocalStokes( double theta, double phi, const cdouble EH[], LocalStokes &locStoke, Stokes &incTransformed );

  /** Saves the Stokes parameters in the phi and theta direction */
  void saveStokesPhiTheta();

  static meep::vec cross( const meep::vec &v1, const meep::vec &v2 );
  static double norm( const meep::vec &vec );

  /** Sets up a rotation matrix around axis 0,1,2 */
  static void setUpRotationMatrix( RotationAxis_t raxis, double alpha, double matrix[3][3] );
  static meep::vec rotateVector( const double rotMat[3][3], const meep::vec &vec );
  static void combineRotationMatrices( const double first[3][3], const double second[3][3], double combined[3][3] );
  static void printRotationMatrix( const double rotMat[3][3] );
};

#endif
