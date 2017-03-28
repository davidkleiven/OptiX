#ifndef VOXELMATERIAL_2D_H
#define VOXELMATERIAL_2D_H
#include <meep.hpp>
#include <armadillo>
#include <string>

struct FileInfo
{
  double voxelsize{1.0};
  unsigned int Nx{0};
  unsigned int Ny{0};
};

struct ComputationalDomain
{
  double xmin{0.0};
  double xmax{1.0};
  double ymin{0.0};
  double ymax{1.0};
};

/** Class for handling input from 2D raw binary files */
class VoxelMaterial2D: public meep::material_function
{
public:
  VoxelMaterial2D(){};

  /** Loads raw voxel file */
  void loadRaw( const char* fname );

  /** Loads raw voxel file */
  void loadRaw( const std::string& fname );

  /** Number of nodes in the x-direction */
  unsigned int Nx() const{ return voxels.n_rows; };

  /** Number of nodes in the y-direction */
  unsigned int Ny() const { return voxels.n_cols; };

  /** Returns the voxel value at index (x,y) */
  unsigned char get( unsigned int x, unsigned int y ) const{ return voxels(x,y); };

  /** Get the voxel value at position r */
  unsigned char get( const meep::vec &r ) const;

  /** True if the vector is inside the structure */
  bool isInside( const meep::vec &r ) const;

  // ====================== PUBLIC ATTRIBUTES ==================================
  /** Threshold when determining if a value is inside the material or not */
  unsigned char threshold{128};

  /** This structure occupies the domain described by domain */
  ComputationalDomain domain;

  /** Holds the information extracted from the filename */
  FileInfo info;
protected:
  static arma::Mat<unsigned char> voxels;

  /** Extracts, resoluion and dimensions from the filename */
  void extractInfoFromFilename( const std::string& fname );

  /** Applies threshold. 1 if the voxel value is larger than threshold, zero otherwise */
  void applyThreshold();

  /** True if the point is outside the computational domain */
  bool isOutsideComputationalDomain( const meep::vec &r ) const;
};

// ================== DISPERSIVE VOXEL 2D MATERIAL =============================
class Voxel2DSusceptibility: public VoxelMaterial2D
{
public:
  Voxel2DSusceptibility( double sigma ): sigma(sigma){};
  Voxel2DSusceptibility( double sigma, double defaultParam ): sigma(sigma), defaultParam(defaultParam){};

  virtual void sigma_row(meep::component c, double sigrow[3], const meep::vec &r) override;
  virtual double chi1p1(meep::field_type ft, const meep::vec &r) override { (void)ft; return f(r); }
  virtual double eps(const meep::vec &r) override { return f(r); }
  virtual double mu(const meep::vec &r) override { return f(r); }
  virtual double conductivity(meep::component c, const meep::vec &r) override {
  (void)c; return f(r); }
  virtual double chi3(meep::component c, const meep::vec &r) override { (void)c; return f(r); }
  virtual double chi2(meep::component c, const meep::vec &r) override { (void)c; return f(r); }
private:
  double f( const meep::vec &r );
  double sigma{0.0};
  double defaultParam{0.0};
};
#endif
