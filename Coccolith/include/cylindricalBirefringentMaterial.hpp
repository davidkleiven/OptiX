#ifndef CYLINDRICAL_BIREFRINGENT_VOXEL_MATERIAL_H
#define CYLINDRICAL_BIREFRINGENT_VOXEL_MATERIAL_H
#include "voxelMaterial.hpp"
#include <armadillo>

struct BirefringentEps
{
  double ordinary{1.0};
  double extraOrdinary{1.0};
};

class CylindricalBirefringentMaterial: public VoxelMaterial
{
public:
  CylindricalBirefringentMaterial(){};

  /** Computes the inertia tensor */
  const arma::mat& inertiaTensor();

  /** Diagonalizes the inertia tensor */
  void diagonalizeInertiaTensor();

  /** Computes the center of mass */
  const arma::vec& centerOfMass();

  /** Returns the rotation axis */
  void rotationAxis( arma::vec &axis );

  /** Returns the principal inertial values */
  const arma::vec& getPrincipalInertia() const { return principalInertia; };

  /** Computes the dielectric tensos. Extraordinary is epsilon along the rotation axis, ordinary is the epsilson perpendicular */
  void dielectricTensor( const meep::vec &r, arma::mat &eps, const BirefringentEps &biref ) const;
private:
  arma::mat inertia;
  arma::vec com;
  arma::vec principalInertia;
  arma::mat principalAxes;
  arma::mat rotationMatrix;
  unsigned int princAxis{0};

  /** Computes the dielectric tensor in the coordinate system spanned by the principal axes */
  void dielectricTensorInPrincCrdSyst( const meep::vec &r, arma::mat &eps, const BirefringentEps &biref ) const;

  /** Transform r in coordinate system to position in the coordinate system spanned by principal axes */
  arma::vec& transformR( const meep::vec &r, arma::vec &transformed ) const;

  /** Sortes the eigenvalues such rotation axis is always the last vector */
  void sortAxes();
};
#endif
