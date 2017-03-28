#ifndef BAND_STRUCTURE_SIMULATION_2D_H
#define BAND_STRUCTURE_SIMULATION_2D_H
#include <string>
#include "voxelMaterial2D.hpp"
#include "refractiveIndexMaterial.hpp"
#include <vector>

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
// ======================== PUBLIC ATTRIBUTES ==================================
  SellmeierMaterial *sellmeier{NULL};
  unsigned int uid{0};
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
  std::vector<double> Ez;

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

};
#endif
