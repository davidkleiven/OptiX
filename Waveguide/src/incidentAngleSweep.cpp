#include "incidentAngleSweep.hpp"

double IncidentAngleSweep::getTheta( unsigned int indx ) const
{
  double step = (theta_max-theta_min)/nTheta;
  return theta_min + step*indx;
}

void IncidentAngleSweep::setWavelength( double wl )
{
  wg.setWavelength( wl );
  pw.setWavelength( wl );
}

void IncidentAngleSweep::setWidth( double width )
{
  wg.setWidth( width );
}

void IncidentAngleSweep::setTransverseDisc( double xmin, double xmax, unsigned int Nx )
{
  double step = (xmax-xmin)/Nx;
  wg.setTransverseDiscretization( xmin, xmax, step );
}

void IncidentAngleSweep::setLongitudinalDisc( double zmin, double zmax, unsigned int Nz )
{
  double step = (zmax-zmin)/Nz;
  wg.setLongitudinalDiscretization( zmin, zmax, step );
}

void IncidentAngleSweep::setCladdingSilicon()
{
  double beta = 1E-7;
  double delta = 1E-5;
  cladding.setRefractiveIndex( delta, beta );
  wg.setCladding( cladding );
}

void IncidentAngleSweep::saveIndx( unsigned int indx )
{
  indxToSave.insert(indx);
}

void IncidentAngleSweep::solve()
{
  solver.setEquation( eq );
  wg.setSolver( solver );
  auto saveIter = indxToSave.begin();
  auto endIter = indxToSave.end();
  for ( unsigned int i=0;i<nTheta;i++ )
  {
    double theta = getTheta( i );
    pw.setAngleDeg(theta);
    wg.setBoundaryConditions( pw );
    wg.solve();
    wg.computeFarField();

    if ( i == 0 )
    {
      farField.setSize( wg.getFarField().n_elem, nTheta );
    }

    for ( unsigned int j=0; j<farField.n_rows;j++ )
    {
      farField(j,i) = wg.getFarField()(j);
    }

    if ( ( saveIter != endIter ) && ( i == *saveIter ) )
    {
      ControlFile ctl("data/incidentAngleSweep");
      wg.save( ctl );
      ctl.save();
      ++saveIter;
    }
  }
}
