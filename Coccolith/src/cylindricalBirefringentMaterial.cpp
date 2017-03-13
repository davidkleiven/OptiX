#include "cylindricalBirefringentMaterial.hpp"
#include <cmath>

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
  rotationMatrix = principalAxes.t();
}

void CylindricalBirefringentMaterial::rotationAxis( arma::vec &axis)
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
    princAxis = 2;
  }
  else if (( abs(I13) < abs(I12) ) && ( abs(I13) < abs(I23) ))
  {
    // 1 and 3 is most similar
    axis = principalAxes.col(1);
    princAxis = 1;
  }
  else
  {
    // 2 and 3 is most similar
    axis = principalAxes.col(0);
    princAxis = 0;
  }
}

void CylindricalBirefringentMaterial::dielectricTensorInPrincCrdSyst( const meep::vec &r, arma::mat &eps, const BirefringentEps &biref ) const
{
  eps.set_size(3,3);
  arma::vec rTrans;
  transformR(r,rTrans);
  for ( unsigned int i=0;i<3;i++ )
  {
    eps(princAxis,i) = 0.0;
    eps(i,princAxis) = 0.0;
  }
  eps(princAxis,princAxis) = biref.ordinary;
  unsigned int indx1 = (princAxis+1)%3;
  unsigned int indx2 = (princAxis+2)%3;
  double R = sqrt( pow(rTrans(indx1),2) + pow(rTrans(indx2),2 ));
  double cosTheta = rTrans(indx1)/R;
  double sinTheta = rTrans(indx2)/R;
  eps(indx1,indx1) = biref.extraOrdinary*cosTheta*cosTheta + biref.ordinary*sinTheta*sinTheta;
  eps(indx2,indx2) = biref.extraOrdinary*sinTheta*sinTheta + biref.ordinary*cosTheta*cosTheta;
  eps(indx1,indx2) = (biref.extraOrdinary - biref.ordinary)*sinTheta*cosTheta;
  eps(indx2,indx1) = eps(indx1,indx2);
}

arma::vec& CylindricalBirefringentMaterial::transformR( const meep::vec &r, arma::vec &transformed ) const
{
  arma::vec original(3);
  original(0) = r.x();
  original(1) = r.y();
  original(2) = r.z();
  transformed = rotationMatrix*original;
  return transformed;
}

void CylindricalBirefringentMaterial::dielectricTensor( const meep::vec &r, arma::mat &eps, const BirefringentEps &biref ) const
{
  dielectricTensorInPrincCrdSyst(r,eps,biref);

  // Rotate the tensor from the axis spanned by the principal values back to the cartesian frame
  eps = principalAxes*eps*rotationMatrix;
}
