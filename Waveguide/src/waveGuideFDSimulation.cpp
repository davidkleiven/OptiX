#include "waveGuideFDSimulation.hpp"
#include "solver2D.hpp"
#include "cladding.hpp"
#include "controlFile.hpp"
#include "crankNicholson.hpp"
#include "waveGuide.hpp"
#include "solver1D.hpp"
#include <H5Cpp.h>
#include <hdf5_hl.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include "paraxialSource.hpp"
#include "borderTracker.hpp"
#include <limits>
//#define DEBUG_BOUNDARY_EXTRACTOR

const double PI = acos(-1.0);

using namespace std;
WaveGuideFDSimulation::WaveGuideFDSimulation(): ParaxialSimulation("noname"){};

WaveGuideFDSimulation::WaveGuideFDSimulation( const char* wgname ): ParaxialSimulation(wgname){}


WaveGuideFDSimulation::~WaveGuideFDSimulation()
{
  if ( wgborder != NULL ) delete wgborder;
  if ( bTracker != NULL ) delete bTracker;
}

void WaveGuideFDSimulation::setCladding( const Cladding &clad )
{
  cladding = &clad;
}

void WaveGuideFDSimulation::setInsideMaterial( const Cladding &clad )
{
  insideMaterial = &clad;
}

void WaveGuideFDSimulation::saveSpecialDatasets( hid_t fid, vector<string> &dsets ) const
{

}

double* WaveGuideFDSimulation::allocateSolutionMatrix() const
{
  double *solution = new double[nodeNumberTransverse()*nodeNumberLongitudinal()];
  return solution;
}

void WaveGuideFDSimulation::deallocateSolutionMatrix( double *matrix ) const
{
  delete [] matrix;
}

void WaveGuideFDSimulation::saveWG( const string &fname ) const
{
  vector<double> xInside;
  vector<double> zInside;

  for ( unsigned int ix=0; ix<nodeNumberTransverse(); ix++ )
  {
    double x = xDisc->min + static_cast<double>(ix)*xDisc->step;
    for ( unsigned int iz=0; iz<nodeNumberLongitudinal(); iz++ )
    {
      double z = zDisc->min + static_cast<double>(iz)*zDisc->step;
      if ( isInsideGuide( x, z ) )
      {
        xInside.push_back(x);
        zInside.push_back(z);
      }
    }
  }

  hsize_t dim = xInside.size();
  unsigned int rank = 1;

  hid_t file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  H5LTmake_dataset( file_id, "xInside", rank, &dim, H5T_NATIVE_DOUBLE, &xInside[0]);
  H5LTmake_dataset( file_id, "zInside", rank, &dim, H5T_NATIVE_DOUBLE, &zInside[0]);
  H5Fclose(file_id);
}

void WaveGuideFDSimulation::extractWGBorders()
{
  wgborder = new vector<WaveGuideBorder>();
  for ( unsigned int iz=0;iz<nodeNumberLongitudinal(); iz++ )
  {
    unsigned int wgNumber = 0;
    bool isInWG = false;
    double z = zDisc->min + iz*zDisc->step;
    if ( waveguideEnded(0.0,z) )
    {
      break;
    }
    for ( unsigned int ix=0;ix<nodeNumberTransverse(); ix++ )
    {
      double x = xDisc->min + ix*xDisc->step;
      if ( isInsideGuide( x, z) )
      {
        if ( !isInWG )
        {
          if ( wgNumber == wgborder->size() )
          {
            WaveGuideBorder border;
            border.x1.push_back(x);
            border.z1.push_back(z);
            wgborder->push_back(border);
          }
          else
          {
            (*wgborder)[wgNumber].x1.push_back(x);
            (*wgborder)[wgNumber].z1.push_back(z);
          }
          isInWG = true;
        }
      }
      else
      {
        if ( isInWG )
        {
          (*wgborder)[wgNumber].x2.push_back(x);
          (*wgborder)[wgNumber].z2.push_back(z);
          isInWG = false;
          wgNumber++;
        }
      }
    }
  }

  // Debugging
  #ifdef DEBUG_BOUNDARY_EXTRACTOR
    for ( auto iter=wgborder->begin(); iter != wgborder->end(); ++iter )
    {
      unsigned int Nx1 = iter->x1.size();
      unsigned int Nz1 = iter->z1.size();
      unsigned int Nx2 = iter->x2.size();
      unsigned int Nz2 = iter->z2.size();
      assert( Nx1 == Nz1 );
      assert( Nx1 == Nx2 );
      assert( Nx1 == Nz2);
    }
  #endif
}

void WaveGuideFDSimulation::sparseSave( const string &fname, double intensityThreshold ) const
{
  vector<double> xSave;
  vector<double> zSave;
  vector<double> intensity;
  unsigned int Nx = nodeNumberTransverse();
  unsigned int Nz = nodeNumberLongitudinal();
  arma::cx_mat sol = solver->getSolution();
  for ( unsigned int ix=0;ix<Nx; ix++ )
  {
    for ( unsigned int iz=0; iz<Nz; iz++ )
    {
      if ( abs(sol(ix, iz)) > intensityThreshold )
      {
        double x = xDisc->min + static_cast<double>(ix)*xDisc->step;
        double z = zDisc->min + static_cast<double>(iz)*zDisc->step;
        xSave.push_back(x);
        zSave.push_back(z);
        intensity.push_back( abs(sol(ix,iz)));
      }
    }
  }

  hsize_t dim = xSave.size();
  hid_t file_id = H5Fcreate( fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  H5LTmake_dataset( file_id, "x", 1, &dim, H5T_NATIVE_DOUBLE, &xSave[0]);
  H5LTmake_dataset( file_id, "z", 1, &dim, H5T_NATIVE_DOUBLE, &zSave[0]);
  H5LTmake_dataset( file_id, "intensity", 1, &dim, H5T_NATIVE_DOUBLE, &intensity[0]);
  H5Fclose(file_id);
}

double WaveGuideFDSimulation::trapezoidalIntegrateIntensityZ( unsigned int iz, unsigned int ixStart, unsigned int ixEnd ) const
{
  double integral = getIntensity( ixStart, iz ) + getIntensity( ixEnd, iz );
  for ( unsigned int ix=ixStart+1; ix <= ixEnd-1; ix ++ )
  {
    double x = getX(ix);
    integral += 2.0*getIntensity( ix, iz );
  }
  double dx = ( getX( ixEnd) - getX( ixStart ) )/static_cast<double>( ixEnd - ixStart );
  return integral*dx*0.5;
}

/*
void WaveGuideFDSimulation::init( const ControlFile &ctl )
{
  // Initialize x-discretization
  xDisc->min = ctl.get()["xDiscretization"]["min"].asDouble();
  xDisc->max = ctl.get()["xDiscretization"]["max"].asDouble();
  xDisc->step = ctl.get()["xDiscretization"]["step"].asDouble();

  // Initialize zDiscretization
  zDisc->min = ctl.get()["zDiscretization"]["min"].asDouble();
  zDisc->max = ctl.get()["zDiscretization"]["max"].asDouble();
  zDisc->step = ctl.get()["zDiscretization"]["step"].asDouble();

  solverInitializedViaInit = true;

  // Need to have a solver object to store the solution.
  // The type should not matter since the solution is already computed.
  // Thus, the solve function should not be called
  solver = new CrankNicholson();
  solver->setGuide( *this );

  // Check that the control file has an entry for the realpart (from early development not all of them do)
  if ( ctl.get().isMember("fieldData") )
  {
    if ( !solver->importHDF5(ctl.get()["datafile"].asString(), ctl.get()["fieldData"].asString()) )
    {
      throw (runtime_error("Error when opening datafile!"));
    }
  }
  else
  {
    // Load the solution
    if ( !solver->importHDF5( ctl.get()["datafile"].asString() ) )
    {
      throw (runtime_error("Error when opening datafile!"));
    }
  }
}
*/

void WaveGuideFDSimulation::getXrayMatProp( double x, double z, double &delta, double &beta) const
{
  // Some assertions for debugging
  assert( cladding != NULL );
  //assert ( x >= xDisc->min );
  //assert ( x <= xDisc->max );
  //assert ( z >= zDisc->min );
  //assert ( z <= zDisc->max );
  double dx = xDisc->step;
  double dz = zDisc->step;

  if ( waveguideEnded(x,z) )
  {
    beta = 0.0;
    delta = 0.0;
    return;
  }

  bool isInside = isInsideGuide( x, z );
  bool neighboursAreInside = isInsideGuide(x+dx,z) && isInsideGuide(x-dx,z) && isInsideGuide(x,z+dz) && \
                             isInsideGuide(x,z-dz);

  double betaInside = 0.0;
  double deltaInside = 0.0;

  if ( insideMaterial != NULL )
  {
    betaInside = insideMaterial->getBeta();
    deltaInside = insideMaterial->getDelta();
  }

  if ( isInside && neighboursAreInside )
  {
    beta = betaInside;
    delta = deltaInside;
    return;
  }
  else if ( isInside && !neighboursAreInside )
  {
    // Average the material properties
    delta = 0.5*( cladding->getDelta() + deltaInside );
    beta =  0.5*( cladding->getBeta() + betaInside );
    return;
  }
  delta = cladding->getDelta();
  beta = cladding->getBeta();
}

cdouble WaveGuideFDSimulation::transverseBC( double z, Boundary_t bnd ) const
{
  assert( cladding != NULL );

  // TODO: Use the getXrayMatProp() function instead
  double delta = cladding->getDelta();
  double beta = cladding->getBeta();
  cdouble im(0.0,1.0);
  double x = 0.0;
  if ( bnd == Boundary_t::TOP )
  {
    x = xDisc->max;
  }
  else if ( bnd == Boundary_t::BOTTOM )
  {
    x = xDisc->min;
  }
  else
  {
    throw (runtime_error("Boundary has to be either TOP or BOTTOM!"));
  }
  return src->get(x,0.0)*exp(-beta*wavenumber*z)*exp(-im*delta*wavenumber*z);
}

cdouble WaveGuideFDSimulation::transverseBC( double z ) const
{
  return this->transverseBC( z, WaveGuideFDSimulation::Boundary_t::BOTTOM );
}

void WaveGuideFDSimulation::useBorderTracker()
{
  if ( bTracker != NULL ) delete bTracker;

  bTracker = new BorderTracker();
  bTracker->setWG(*this);
  bTracker->init();
}

void WaveGuideFDSimulation::save( ControlFile &ctl )
{
  ParaxialSimulation::save( ctl );
  if ( wgborder != NULL )
  {
    unsigned int indx = 0;
    for ( auto iter=wgborder->begin(); iter != wgborder->end(); ++iter )
    {
      stringstream uxname, uzname, bxname, bzname;
      uxname << "upperBorderX" << indx;
      uzname << "upperBorderZ" << indx;
      bxname << "lowerBorderX" << indx;
      bzname << "lowerBorderZ" << indx;
      indx++;
      hsize_t dim = iter->x1.size();
      const double *ptr = &iter->x1[0];
      arma::vec vec(iter->x1);
      unsigned int rank = 1;
      dsetnames.push_back( bxname.str() );
      saveVec( iter->x1, dsetnames.back().c_str() );
      dsetnames.push_back( bzname.str() );
      saveVec( iter->z1, dsetnames.back().c_str() );
      dsetnames.push_back( uxname.str() );
      saveVec( iter->x2, dsetnames.back().c_str() );
      dsetnames.push_back( uzname.str() );
      saveVec( iter->z2, dsetnames.back().c_str() );
    }
  }
}

cdouble WaveGuideFDSimulation::padExitField( double x, double z ) const
{
  double delta = cladding->getDelta();
  double beta = cladding->getBeta();
  cdouble im(0.0,1.0);
  return src->get(x,0.0)*exp(-beta*wavenumber*z)*exp(-im*delta*wavenumber*z);
}
