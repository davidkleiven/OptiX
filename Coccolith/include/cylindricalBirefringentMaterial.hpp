#ifndef CYLINDRICAL_BIREFRINGENT_VOXEL_MATERIAL_H
#define CYLINDRICAL_BIREFRINGENT_VOXEL_MATERIAL_H
#include "voxelMaterial.hpp"
#include <armadillo>

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
  void rotationAxis( arma::vec &axis ) const;

  /** Returns the principal inertial values */
  const arma::vec& getPrincipalInertia() const { return principalInertia; };
private:
  arma::mat inertia;
  arma::vec com;
  arma::vec principalInertia;
  arma::mat principalAxes;
};
#endif
