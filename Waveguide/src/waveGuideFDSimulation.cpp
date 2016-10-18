#include "waveGuideFDSimulation.hpp"
#include "solver2D.hpp"
#include "cladding.hpp"
#include "controlFile.hpp"
#include "crankNicholson.hpp"
#include <H5Cpp.h>
#include <hdf5_hl.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>

const double PI = acos(-1.0);

using namespace std;
WaveGuideFDSimulation::WaveGuideFDSimulation(): xDisc(new Disctretization), zDisc(new Disctretization), name(""){};

WaveGuideFDSimulation::WaveGuideFDSimulation( const char* wgname ): WaveGuideFDSimulation()
{
  name = wgname;
}

void WaveGuideFDSimulation::solve()
{
  if ( solver == NULL )
  {
    throw ( runtime_error("No solver specified!"));
  }
  solver->solve();
}

WaveGuideFDSimulation::~WaveGuideFDSimulation()
{
  delete xDisc;
  delete zDisc;

  if ( solverInitializedViaInit )
  {
    delete solver;
  }
}

void WaveGuideFDSimulation::setWaveLength( double lambda )
{
  wavenumber = 2.0*PI/lambda;
}

void WaveGuideFDSimulation::setTransverseDiscretization( double xmin, double xmax, double step )
{
  xDisc->min = xmin;
  xDisc->max = xmax;
  xDisc->step = step;
}

void WaveGuideFDSimulation::setLongitudinalDiscretization( double xmin, double xmax, double step )
{
  zDisc->min = xmin;
  zDisc->max = xmax;
  zDisc->step = step;
}

unsigned int WaveGuideFDSimulation::nodeNumberTransverse() const
{
  return (xDisc->max - xDisc->min)/xDisc->step + 1.0;
}

unsigned int WaveGuideFDSimulation::nodeNumberLongitudinal() const
{
  return ( zDisc->max - zDisc->min)/zDisc->step + 1.0;
}

void WaveGuideFDSimulation::setSolver( Solver2D &solv )
{
  if ( solverInitializedViaInit )
  {
    throw ("You cannot set a new solver when it has been initialized via the init function");
  }

  solver = &solv;
  solver->setGuide( *this );
}

void WaveGuideFDSimulation::setCladding( const Cladding &clad )
{
  cladding = &clad;
}

void WaveGuideFDSimulation::save( ControlFile &ctl ) const
{
  save(ctl, -1.0);
}

void WaveGuideFDSimulation::save( ControlFile &ctl, double intensityThreshold ) const
{
  if ( solver == NULL )
  {
    throw ( runtime_error("No solver specified!\n"));
  }

  bool useSparse = (intensityThreshold > 0.0);
  string fname = ctl.getFnameTemplate();
  string h5fname = fname+".h5";
  string h5fieldfname = fname+"_field.h5";
  string jsonfname = fname+".json";
  string wgFname = fname+"_wg.h5";

  //arma::abs(solver->getSolution()).save(h5fname.c_str(), arma::hdf5_binary);
  if ( useSparse )
  {
    sparseSave( h5fname, intensityThreshold );
    clog << "Solution written to " << h5fname << endl;
  }
  else
  {
    arma::mat absSol = arma::abs(solver->getSolution());
    absSol.save(h5fname.c_str(), arma::hdf5_binary);
    clog << "Amplitude written to " << h5fname << endl;

    absSol = arma::real(solver->getSolution());
    absSol.save(h5fieldfname.c_str(), arma::hdf5_binary);
    clog << "Realpart written to " << h5fieldfname << endl;
  }

  saveWG( wgFname );
  clog << "Points inside waveguide written to " << wgFname << endl;

  Json::Value wginfo;
  Json::Value solverInfo;
  wginfo["Cladding"]["delta"] = cladding->getDelta();
  wginfo["Cladding"]["beta"] = cladding->getBeta();
  ctl.get()["datafile"] = h5fname;
  ctl.get()["fieldData"] = h5fieldfname;
  ctl.get()["wgfile"] = wgFname;
  ctl.get()["name"] = name;
  ctl.get()["sparseSave"] = useSparse;
  ctl.get()["sparseThreshold"] = intensityThreshold;
  ctl.get()["xDiscretization"]["min"] = xDisc->min;
  ctl.get()["xDiscretization"]["max"] = xDisc->max;
  ctl.get()["xDiscretization"]["step"] = xDisc->step;
  ctl.get()["zDiscretization"]["min"] = zDisc->min;
  ctl.get()["zDiscretization"]["max"] = zDisc->max;
  ctl.get()["zDiscretization"]["step"] = zDisc->step;
  fillInfo( wginfo );
  // TODO: For some reason the next line gives a segmentation fault
  //solver->fillInfo( solverInfo );
  ctl.get()["solver"] = solverInfo;
  ctl.get()["waveguide"] = wginfo;
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

void WaveGuideFDSimulation::closestIndex( double x, double z, unsigned int &ix, unsigned int &iz ) const
{
  ix = (x - xDisc->min)/xDisc->step;
  iz = (z - zDisc->min)/zDisc->step;
}

double WaveGuideFDSimulation::getIntensity( double x, double z ) const
{
  // Some assertions for debugging
  assert( x < xDisc->max );
  assert( x >= xDisc->min );
  assert( z < zDisc->max );
  assert( z >= zDisc->min );

  unsigned int ix, iz;
  closestIndex( x, z, ix, iz );

  double x1 = xDisc->min + ix*xDisc->step;
  double x2 = x1 + xDisc->step;
  double z1 = zDisc->min + iz*zDisc->step;
  double z2 = z1 + iz*zDisc->step;
  arma::cx_mat sol = solver->getSolution();

  double intensity = pow( abs(sol(ix,iz)),2 )*(x2 - x)*(z2-z);
  intensity += pow( abs(sol(ix+1,iz)), 2 )*(x-x1)*(z2-z);
  intensity += pow( abs(sol(ix,iz+1)), 2 )*(z-z1)*(x2-x);
  intensity += pow( abs(sol(ix+1,iz+1)), 2 )*(z2-z)*(x-x1);
  return intensity/( (x2-x1)*(z2-z1) );
}

double WaveGuideFDSimulation::getIntensity( unsigned int ix, unsigned int iz ) const
{
  return pow( abs( solver->getSolution()(ix,iz) ), 2 );
}

double WaveGuideFDSimulation::getZ( unsigned int iz ) const
{
  return zDisc->min + iz*zDisc->step;
}

double WaveGuideFDSimulation::getX( unsigned int ix ) const
{
  return xDisc->min + ix*xDisc->step;
}

double WaveGuideFDSimulation::trapezoidalIntegrateIntensityZ( unsigned int iz, unsigned int ixStart, unsigned int ixEnd ) const
{
  double integral = getIntensity( ixStart, iz ) + getIntensity( ixEnd, iz );
  for ( unsigned int ix=ixStart+1; ix <= ixEnd-1; ix ++ )
  {
    integral += 2.0*getIntensity( ix, iz );
  }
  double dx = ( getX( ixEnd) - getX( ixStart ) )/static_cast<double>( ixEnd - ixStart );
  return integral*dx*0.5;
}

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
