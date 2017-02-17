#ifndef VOXEL_MATERIAL_H
#define VOXEL_MATERIAL_H
#include <meep.hpp>
#include <armadillo>
#include <string>

/** Read info from filename */
struct InfoFromFilename
{
  double voxelsize{0.0};
  unsigned int Nx{0};
  unsigned int Ny{0};
  unsigned int Nz{0};
};

class VoxelMaterial: public meep::material_function
{
public:
  VoxelMaterial(){};

  /** Returns the dielectric function as a function of position */
  virtual double eps( const meep::vec &r ) override;

  /** Returns the conductivity as a function of position */
  virtual double conductivity( meep::component c, const meep::vec &r ) override;

  /** Loads data from raw file without header */
  void loadRaw( const char* fname );

  /** Loads data from raw file */
  void loadRaw( const std::string &fname );
private:
  /** Extracts the dimension of the data from the filename */
  static void extractDimsFromFilename( const std::string &fname, InfoFromFilename &info );

  /** Separates the data in two sets */
  void applyThreshold();

  /** Print statistics of voxels */
  void showStatistics() const;

  arma::Cube<unsigned char> voxels;
};

#endif
