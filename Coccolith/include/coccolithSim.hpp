#ifndef COCCOLITH_SIMULATION_H
#define COCCOLITH_SIMULATION_H
#include <string>
#include "voxelMaterial.hpp"
#include "fieldMonitors.hpp"
#include <complex>
#include <visa/visa.hpp>
#include <cstdlib>
#define UID_MAX 10000000

typedef std::complex<double> cdouble;
enum class MainPropDirection_t{ X, Y, Z };
typedef MainPropDirection_t IntegrationDir_t;
enum class Plane_t{XY, XZ, YZ};

enum class SourcePosition_t{TOP, BOTTOM};

class CoccolithSimulation
{
public:
  CoccolithSimulation(){};
  virtual ~CoccolithSimulation();

  /** Load voxel model from raw binary file */
  void loadVoxels( const char* fname );

  /** Set the incident wavevector. Note that this is a static member */
  void setIncWaveVector( const meep::vec &wave ){ waveVec = wave; };

  /** Set the main propagation direction */
  void setMainPropagationDirection( MainPropDirection_t propDir );

  /** Return the wavelength in pixel units */
  double getWavelength() const;

  /** Set the PML thickness in number of wavelengths */
  void setPMLInWavelengths( double newThick );

  /** Run the simulation */
  void run();

  /** Initialize the source profile */
  void initSource( double freq, double fwidth );

  /** Initialize the simulation */
  void init();

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
  void exportResults();

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
private:
  CaCO3Cocco material;
  MainPropDirection_t propagationDir{MainPropDirection_t::Z};

  meep::volume* srcVol{NULL};
  meep::structure* struc{NULL};
  meep::grid_volume gdvol;
  meep::fields* field{NULL};
  meep::src_time *sourceTime{NULL}; // Do not delete this pointer
  meep::gaussian_src_time *source{NULL};
  unsigned int uid{0};
  std::string outdir{"data/"};
  unsigned int nSave{30};
  bool isInitialized{false};
  unsigned int plotUpdateFreq{30};

  double pmlThicknessInWavelengths{3.0};
  double centerFrequency{1.0};
  double freqWidth{0.5};
  double resolution{1.0};
  bool materialLoaded{false};
  unsigned int nfreq{100};
  unsigned int nMonitorX{256};
  unsigned int nMonitorY{256};
  unsigned int nMonitorZ{256};
  double tEnd{100.0};
  visa::WindowHandler *plots{NULL};
  bool userOverridedEndTime{false};
  bool realTimeVisualization{true};
  bool geoIsInitialized{false};

  /** Visualized intensity */
  void visualize();

  // Flux planes
  meep::volume *dftVolTransmit{NULL};
  meep::dft_flux *transmitFlux{NULL};

  // Monitor planes 1
  FieldMonitor *monitor1{NULL};
  FieldMonitor *monitor2{NULL};

  arma::mat bkg1;
  arma::mat bkg2;


  // Source corners
  meep::vec crn1;
  meep::vec crn2;

  SourcePosition_t srcPos{SourcePosition_t::TOP};

  meep::h5file *file{NULL};

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
};

#endif
