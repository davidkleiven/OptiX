#include <cstdlib>
#include <libscuff.h>

const double PI = acos(-1.0);

void visualize( double ymin, double ymax, double zmin, double zmax, unsigned int nPoints, scuff::RWGGeometry &geo, \
                double omega, double kBloch[2], IncField &IF, HVector &rhs )
{

}

       
  
int main(int argc, char **argv)
{
  scuff::RWGGeometry::AssignBasisFunctionsToExteriorEdges=false;
  scuff::RWGGeometry geo = scuff::RWGGeometry("halfspace.scuffgeo");
  SetLogFileName("fresnel.log");
  geo.SetLogLevel(SCUFF_VERBOSELOGGING);
 
  // Allocate BEM matrix and rhs vector
  HMatrix *matrix = geo.AllocateBEMMatrix();
  HVector *rhsVec = geo.AllocateRHSVector();

  // Source definition
  double angle = 40.0*PI/180.0;
  double kHat[3] = {0.0,sin(angle),cos(angle)};
  cdouble E0_s[3]={1.0,0.0,0.0};
  PlaneWave pw(E0_s, kHat);
  double omega = 1.0;
  
  double kBloch[2] = {0.0,0.0};
 
  // Assembling BEM matrix
  double ksource = real(geo.RegionMPs[pw.RegionIndex]->GetRefractiveIndex(omega))*omega;
  kBloch[1] = ksource*sin(angle);

  geo.AssembleBEMMatrix(static_cast<cdouble>(omega), kBloch, matrix);
  matrix->LUFactorize();

  // Solve for s polarisation
  geo.AssembleRHSVector(static_cast<cdouble>(omega), kBloch, &pw, rhsVec);
  int info = matrix->LUSolve(rhsVec);

  // Fill evaluation points
  unsigned int nPoints = 4;
  HMatrix Xpoints(nPoints*nPoints, 3);
  double z = 0.0;
  double dz = 1.0/static_cast<double>( nPoints );
  double dy = 1.0/static_cast<double>( nPoints );
  for ( unsigned int iz=0;iz<nPoints;iz++)
  {
    double y = 0.0;
    for ( unsigned int iy=0;iy<nPoints;iy++)
    {
      Xpoints.SetEntry( iz*nPoints+iy, 0, 0.5);
      Xpoints.SetEntry( iz*nPoints+iy, 1, y );
      Xpoints.SetEntry( iz*nPoints+iy, 2, z );
      y += dy;
    }
    z += dz;
  }
  HMatrix field( nPoints*nPoints, 6, LHM_COMPLEX );
  
  // Get fields
  // 1) Only incident
  //geo.GetFields( &pw, NULL, omega, kBloch, &Xpoints, &field );

  // 2) Only scattered
  geo.GetFields( NULL, rhsVec, omega, kBloch, &Xpoints, &field );

  // 3) Total
  //geo.GetFields( &pw, rhsVec, omega, kBloch, &Xpoints, &field );
  delete matrix;
  delete rhsVec;
  return 0;
}
