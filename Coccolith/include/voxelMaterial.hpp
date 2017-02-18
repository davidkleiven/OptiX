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

  /** Loads data from raw file without header */
  void loadRaw( const char* fname );

  /** Loads data from raw file */
  void loadRaw( const std::string &fname );

  /** Shows an animation where it slides through the voxel array in different direction */
  void slideThroughVoxelArray() const;

  /** Shows projections in all three directions */
  void showProjections() const;

  /** Domain size in the X-direction */
  unsigned int sizeX() const { return voxels.n_rows; };

  /** Domain size in the Y-direction */
  unsigned int sizeY() const { return voxels.n_cols; };

  /** Domain size in the Z-direction */
  unsigned int sizeZ() const { return voxels.n_slices; };
protected:
  /** Extracts the dimension of the data from the filename */
  static void extractDimsFromFilename( const std::string &fname, InfoFromFilename &info );

  static void fillArmaMat( const arma::Mat<unsigned char> &values, arma::mat &matrix );

  /** Separates the data in two sets */
  void applyThreshold();

  /** Print statistics of voxels */
  void showStatistics() const;

  arma::Cube<unsigned char> voxels;
};

class CaCO3Cocco: public VoxelMaterial
{
public:
  CaCO3Cocco(){};

  /** Returns the dielectric function as a function of position */
  virtual double eps( const meep::vec &r ) override;

  /** Returns the conductivity as a function of position */
  virtual double conductivity( meep::component c, const meep::vec &r ) override;

private:
  double epsilon{2.72};
};

#endif
