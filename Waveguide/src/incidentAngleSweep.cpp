#include "incidentAngleSweep.hpp"
#include "controlFile.hpp"
#include <iostream>
#include <H5Cpp.h>
#include <hdf5_hl.h>
#include <cstdlib>
#include <ctime>
#include <sstream>

using namespace std;
IncidentAngleSweep::IncidentAngleSweep():picDir("")
{
  srand( time(0) );
  uid = rand()%1000000;
}

double IncidentAngleSweep::getAngle( unsigned int indx ) const
{
  double step = (thetaMax-thetaMin)/nTheta;
  return thetaMin + step*indx;
}

void IncidentAngleSweep::setWavelength( double wl )
{
  wg.setWaveLength( wl );
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

void IncidentAngleSweep::setIncAngles( double min, double max, unsigned int N )
{
  thetaMin = min;
  thetaMax = max;
  nTheta = N;
}
void IncidentAngleSweep::solve()
{
  solver.setEquation( eq );
  wg.setSolver( solver );
  auto saveIter = indxToSave.begin();
  auto endIter = indxToSave.end();
  for ( unsigned int i=0;i<nTheta;i++ )
  {
    clog << "Running "<< i+1 << " of " << nTheta << endl;
    double theta = getAngle( i );
    pw.setAngleDeg(theta);
    wg.setBoundaryConditions( pw );
    wg.solve();
    wg.computeFarField(fftSignalLength);

    if ( i == 0 )
    {
      farField.set_size( wg.getFarField().n_elem, nTheta );
    }

    for ( unsigned int j=0; j<farField.n_rows;j++ )
    {
      farField(j,i) = wg.getFarField()(j);
    }

    if ( ( saveIter != endIter ) && ( i == *saveIter ) )
    {
      ControlFile ctl("data/incidentAngleSweep");
      wg.extractWGBorders();
      wg.save( ctl );
      ctl.save();
      ++saveIter;
    }

    if ( vis != NULL )
    {
      vis->fillVertexArray( arma::abs( wg.getSolver().getSolution() ) );
      vis->display();

      if ( picDir != "" )
      {
        auto image = vis->capture();
        stringstream ss;
        ss << picDir << "/img" << i << ".png";
        image.saveToFile(ss.str().c_str());
      }
    }
  }
  vis->close();
}

void IncidentAngleSweep::save( const string &fname ) const
{
    stringstream ss;
    ss << fname << uid << ".h5";
    const double PI = acos(-1.0);
    double qMax = PI/wg.transverseDiscretization().step;
    double phiMax = qMax/wg.getWavenumber();
    phiMax *= 180.0/PI;
    double phiMin = -phiMax;

    hid_t file_id = H5Fcreate( ss.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dim[2] = {farField.n_cols, farField.n_rows};
    H5LTmake_dataset( file_id, "intensity", 2, dim, H5T_NATIVE_DOUBLE, farField.memptr());
    H5LTset_attribute_double( file_id, "intensity", "thetaMin", &thetaMin, 1);
    H5LTset_attribute_double( file_id, "intensity", "thetaMax", &thetaMax, 1);
    H5LTset_attribute_double( file_id, "intensity", "phiMin", &phiMin, 1);
    H5LTset_attribute_double( file_id, "intensity", "phiMax", &phiMax, 1);
    H5Fclose( file_id );
    clog << "Results written to " << ss.str() << endl;
  }

void IncidentAngleSweep::setAlcoholInside()
{
  double delta = 4.08E-6; // For 7.9 keV (0.1569 nm)
  inside.setRefractiveIndex(delta, 0.0);
  wg.setInsideMaterial(inside);
}

void IncidentAngleSweep::savePic( const char *dir )
{
  picDir = dir;
}
