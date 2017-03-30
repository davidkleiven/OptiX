#ifndef BAND_STRUCTURE_SIMULATION_2D_H
#define BAND_STRUCTURE_SIMULATION_2D_H
#include <string>
#include "voxelMaterial2D.hpp"
#include "refractiveIndexMaterial.hpp"
#include <vector>
#include <complex>
#include <armadillo>
typedef std::complex<double> cdouble;

class HarminvRes
{
public:
  HarminvRes( unsigned int nfreq );
  arma::Col<cdouble> amplitude;
  arma::vec freqRe;
  arma::vec freqIm;
  arma::vec freqErr;
  int numModesFound{0};
};

class BandStructure2DSimulation
{
public:
  BandStructure2DSimulation():bloch(0.0,0.0){};
  ~BandStructure2DSimulation();

  /** Sets the material and updates initializes gdvol */
  void setMaterial( VoxelMaterial2D &material );

  /** Runs the simulation */
  void run();

  /** Save the results */
  void save();

  /** Analyse the results with harminv */
  void findModes();
// ======================== PUBLIC ATTRIBUTES ==================================
  SellmeierMaterial *sellmeier{NULL};
  unsigned int uid{0};
  unsigned int nHarminvFreq{100};
  double resolution{1.0};
  double freq{1.0};
  double freqWidth{0.5};
  unsigned int nfreq{200};
  meep::vec bloch;
  bool useSingleFrequency{false};
  bool addTimestampToFilename{true};
  bool addUIDToFilename{false};
  std::string prefix{"bandStructureDefault"};
private:
  std::string timestamp;
  meep::grid_volume gdvol;
  meep::src_time *srcTime{NULL};
  meep::dft_ldos *ldos{NULL};
  meep::vec srcPos;
  std::vector<cdouble> Ez;

  VoxelMaterial2D *material{NULL};
  meep::h5file *outfile{NULL};

  /** Initialize the simulation */
  void init();

  /** Updates the structure with the parameters in the Sellmeier material */
  void updateStructure();

  /** Prints info about the simulation */
  void printInfo() const;

  /** Save parameters */
  void saveParams();

  /** Initializes the timestamp string */
  void setTimeStamp();

  meep::structure *struc{NULL};
  meep::fields *field{NULL};
  HarminvRes *modes{NULL};

};
#endif
