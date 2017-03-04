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

struct DomainSize
{
  double xmin, xmax;
  double ymin, ymax;
  double zmin, zmax;
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

  /** Get profile projection onto the XY-plane */
  void projectionXY( arma::mat &mat ) const;

  /** Get the profile projected onto the XZ-plane */
  void projectionXZ( arma::mat &mat ) const;

  /** Get the profile projected onto the YZ-plane */
  void projectionYZ( arma::mat &mat ) const;

  /** Domain size in the X-direction */
  unsigned int sizeX() const { return voxels.n_rows; };

  /** Domain size in the Y-direction */
  unsigned int sizeY() const { return voxels.n_cols; };

  /** Domain size in the Z-direction */
  unsigned int sizeZ() const { return voxels.n_slices; };

  /** Returns the voxel size in nano meters */
  double getVoxelSize() const{ return vxsize; };

  /** Set the domain size. NOTE: Assumed to correspond to the array */
  void setDomainSize( const meep::grid_volume &gvol, double PMLThickInPx );

  /** If set to true the refractive index is 1 in the entire domain */
  void setReferenceRun( bool newval ){ referenceRun = newval; };

  /** True if the run is a reference run */
  bool isReferenceRun() const{ return referenceRun; };
protected:
  /** Extracts the dimension of the data from the filename */
  static void extractDimsFromFilename( const std::string &fname, InfoFromFilename &info );

  static void fillArmaMat( const arma::Mat<unsigned char> &values, arma::mat &matrix );

  /** Separates the data in two sets */
  void applyThreshold();

  /** Print statistics of voxels */
  void showStatistics() const;

  /** Translate meep coordinates to index in voxel array */
  void meepVecToIndx( const meep::vec &r, unsigned int indx[3] );

  DomainSize domain;

  double vxsize{1.0};

  arma::Cube<unsigned char> voxels;

  bool referenceRun{false};

  /** Returns true if the position inside the domain */
  bool isInsideDomain( const meep::vec &r ) const;
};

class CaCO3Cocco: public VoxelMaterial
{
public:
  CaCO3Cocco(){};

  /** Returns the dielectric function as a function of position */
  virtual double eps( const meep::vec &r ) override;

  /** Returns the conductivity as a function of position */
  virtual double conductivity( meep::component c, const meep::vec &r ) override;

  /** MEEP function that has to implemented in order to get it to work */
  virtual double chi1p1(meep::field_type ft, const meep::vec &r) { return eps(r); }

private:
  double epsilon{2.72};
};

#endif
