#include "cylindricalBirefringentMaterial.hpp"

const arma::vec& CylindricalBirefringentMaterial::centerOfMass()
{
  com.set_size(3);
  com.fill(0.0);
  double mass = 0.0;
  for ( unsigned int z=0;z<Nz();z++ )
  for ( unsigned int y=0;y<Ny();y++ )
  for ( unsigned int x=0;x<Nx();x++ )
  {
    com(0) += x*get(x,y,z);
    com(1) += y*get(x,y,z);
    com(2) += z*get(x,y,z);
    mass += get(x,y,z);
  }
  com /= mass;
  return com;
}

const arma::mat& CylindricalBirefringentMaterial::inertiaTensor()
{
  centerOfMass();
  inertia.set_size(3,3);
  inertia.fill(0.0);
  for ( unsigned int z=0;z<Nz();z++ )
  for ( unsigned int y=0;y<Ny();y++ )
  for ( unsigned int x=0;x<Nx();x++ )
  {
    double xcrd = static_cast<int>(x)-com(0);
    double ycrd = static_cast<int>(y)-com(1);
    double zcrd = static_cast<int>(z)-com(2);
    inertia(0,0) += (ycrd*ycrd + zcrd*zcrd)*get(x,y,z);
    inertia(0,1) -= xcrd*ycrd*get(x,y,z);
    inertia(0,2) -= xcrd*zcrd*get(x,y,z);
    inertia(1,1) += (xcrd*xcrd+zcrd*zcrd)*get(x,y,z);
    inertia(1,2) -= ycrd*zcrd*get(x,y,z);
    inertia(2,2) += (xcrd*xcrd+ycrd*ycrd)*get(x,y,z);
  }
  inertia(1,0) = inertia(0,1);
  inertia(2,1) = inertia(1,2);
  inertia(2,0) = inertia(0,2);
  return inertia;
}

void CylindricalBirefringentMaterial::diagonalizeInertiaTensor()
{
  arma::eig_sym( principalInertia, principalAxes, inertia );
}

void CylindricalBirefringentMaterial::rotationAxis( arma::vec &axis) const
{
  // If there is some kind of cylindrical symmetry, two of the principal
  // inertia should be equal.
  double I12 = principalInertia(0)-principalInertia(1);
  double I13 = principalInertia(0)-principalInertia(2);
  double I23 = principalInertia(1)-principalInertia(2);

  if (( abs(I12) < abs(I13) ) && ( abs(I12) < abs(I23) ))
  {
    // 1 and 2 is most similar
    axis = principalAxes.col(2);
  }
  else if (( abs(I13) < abs(I12) ) && ( abs(I13) < abs(I23) ))
  {
    // 1 and 3 is most similar
    axis = principalAxes.col(1);
  }
  else
  {
    // 2 and 3 is most similar
    axis = principalAxes.col(0);
  }
}
