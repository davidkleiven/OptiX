#include "cylindricalBirefringentMaterial.hpp"
#include <iostream>
#include <string>
#include <armadillo>
#include <mpi.h>

using namespace std;

int main( int argc, char** argv )
{
  meep::initialize mpi(argc, argv);
  CylindricalBirefringentMaterial material;
  material.loadRaw( "data/cocco8cv4Rotated_216_182_249_253.raw" );

  cout << "COM:\n" << material.centerOfMass() << endl;
  cout << "Inertia tensor:\n" << material.inertiaTensor() << endl;
  material.diagonalizeInertiaTensor();
  cout << "Principal inertia:\n" << material.getPrincipalInertia() << endl;
  arma::vec axis;
  material.rotationAxis(axis);
  cout << "Rotation axis:\n" << axis << endl;

  BirefringentEps biref;
  biref.ordinary = 1.63;
  biref.extraOrdinary = 1.3;

  arma::mat epsilon;
  meep::vec r(1.0,0.0,0.0);
  material.dielectricTensor( r, epsilon, biref );
  cout << "Dielectric tensor:\n" << epsilon << endl;
  return 0;
}
