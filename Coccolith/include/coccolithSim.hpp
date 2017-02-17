#ifndef COCCOLITH_SIMULATION_H
#define COCCOLITH_SIMULATION_H
#include <string>
#include "voxelMaterial.hpp"
#include <complex>

typedef std::complex<double> cdouble;
enum class MainPropDirection_t{ X, Y, Z };

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

  /** Initialize the source profile */
  void initSource( double freq, double fwidth );

  /** Initialize the simulation */
  void init();

  /** Set the number of frequencies to be used in the harmonic inversion */
  void setNfreqFT( unsigned int numberOfFreq ){ nfreq = numberOfFreq; };

  /** Set the number of timesteps to save */
  void setNumberOfImgSave( unsigned int newnSave){ nSave = newnSave; };
private:
  CaCO3Cocco material;
  MainPropDirection_t propagationDir{MainPropDirection_t::Z};

  meep::volume* srcVol{NULL};
  meep::structure* struc{NULL};
  meep::grid_volume gdvol;
  meep::fields* field{NULL};
  meep::gaussian_src_time *source{NULL};
  unsigned int uid{0};
  std::string outdir{"data/"};
  unsigned int nSave{30};

  double pmlThicknessInWavelengths{3.0};
  double centerFrequency{1.0};
  double freqWidth{0.5};
  double resolution{1.0};
  bool materialLoaded{false};
  unsigned int nfreq{100};

  // Flux planes
  meep::volume *dftVolSource{NULL};
  meep::volume *dftVolTransmit{NULL};
  meep::dft_flux *srcFlux{NULL};
  meep::dft_flux *transmitFlux{NULL};

  /** Add a source volume */
  void addSourceVolume();

  /** Add a source to the field */
  void addSource();

  /** Add a structure for the simulation */
  void addStructure();

  /** Add fields */
  void addFields();

  /** Add DFT flux planes */
  void addFluxPlanes( const meep::vec &srcCrn1, const meep::vec &srcCrn2 );

  /** Returns the PML thickness in pixel units */
  double getPMLThickness() const;

  static meep::vec waveVec;

  static cdouble amplitude( const meep::vec &r );
};

#endif
