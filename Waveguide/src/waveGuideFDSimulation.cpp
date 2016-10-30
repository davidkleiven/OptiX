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
//#define DEBUG_BOUNDARY_EXTRACTOR

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

  if ( src == NULL )
  {
    throw ( runtime_error("No source specified!"));
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

  if ( farFieldModulus != NULL ) delete farFieldModulus;
  if ( wgborder != NULL ) delete wgborder;
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

  if ( src == NULL )
  {
    throw ( runtime_error("No source specified!\n") );
  }
  bool useSparse = (intensityThreshold > 0.0);
  string fname = ctl.getFnameTemplate();
  string h5fname = fname+".h5";
  string h5fieldfname = fname+"_field.h5";
  string jsonfname = fname+".json";
  string wgFname = fname+"_wg.h5";
  string farFieldFname = fname+"_farField.h5";
  string phaseName = fname+"_phase.h5";

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

    solver->getField( absSol );
    absSol.save(h5fieldfname.c_str(), arma::hdf5_binary);
    clog << "Realpart written to " << h5fieldfname << endl;

    solver->getPhase( absSol );
    absSol.save(phaseName.c_str(), arma::hdf5_binary);
    clog << "Phase written to " << phaseName << endl;
  }

  hid_t file_id = H5Fcreate(wgFname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
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
      unsigned int rank = 1;
      H5LTmake_dataset( file_id, bxname.str().c_str(), rank, &dim, H5T_NATIVE_DOUBLE, ptr);
      ptr = &iter->z1[0];
      H5LTmake_dataset( file_id, bzname.str().c_str(), rank, &dim, H5T_NATIVE_DOUBLE, ptr);
      ptr = &iter->x2[0];
      H5LTmake_dataset( file_id, uxname.str().c_str(), rank, &dim, H5T_NATIVE_DOUBLE, ptr);
      ptr = &iter->z2[0];
      H5LTmake_dataset( file_id, uzname.str().c_str(), rank, &dim, H5T_NATIVE_DOUBLE, ptr);
    }
    H5Fclose(file_id);
    ctl.get()["wgfile"] = wgFname;
    clog << "Waveguide borders written to " << wgFname << endl;
  }

  Json::Value wginfo;
  Json::Value solverInfo;
  Json::Value sourceInfo;
  src->info( sourceInfo );
  wginfo["Cladding"]["delta"] = cladding->getDelta();
  wginfo["Cladding"]["beta"] = cladding->getBeta();
  wginfo["length"] = wglength;
  ctl.get()["datafile"] = h5fname;
  ctl.get()["fieldData"] = h5fieldfname;
  ctl.get()["phase"] = phaseName;
  ctl.get()["name"] = name;
  ctl.get()["sparseSave"] = useSparse;
  ctl.get()["sparseThreshold"] = intensityThreshold;
  ctl.get()["xDiscretization"]["min"] = xDisc->min;
  ctl.get()["xDiscretization"]["max"] = xDisc->max;
  ctl.get()["xDiscretization"]["step"] = xDisc->step;
  ctl.get()["zDiscretization"]["min"] = zDisc->min;
  ctl.get()["zDiscretization"]["max"] = zDisc->max;
  ctl.get()["zDiscretization"]["step"] = zDisc->step;
  ctl.get()["source"] = sourceInfo;
  fillInfo( wginfo );
  // TODO: For some reason the next line gives a segmentation fault
  //solver->fillInfo( solverInfo );
  ctl.get()["solver"] = solverInfo;
  ctl.get()["waveguide"] = wginfo;

  if ( farFieldModulus != NULL )
  {
    ctl.get()["farFieldFile"] = farFieldFname;
    saveFarField( farFieldFname, ctl.getUID() );
  }
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
    double x = getX(ix);
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

void WaveGuideFDSimulation::getXrayMatProp( double x, double z, double &delta, double &beta) const
{
  // Some assertions for debugging
  assert( cladding != NULL );
  assert ( x >= xDisc->min );
  assert ( x <= xDisc->max );
  assert ( z >= zDisc->min );
  assert ( z <= zDisc->max );
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
  if ( isInside && neighboursAreInside )
  {
    beta = 0.0;
    delta = 0.0;
    return;
  }
  else if ( isInside && !neighboursAreInside )
  {
    // Average the material properties
    delta = 0.5*cladding->getDelta();
    beta =  0.5*cladding->getBeta();
    return;
  }
  delta = cladding->getDelta();
  beta = cladding->getBeta();
}

void WaveGuideFDSimulation::getExitField( arma::vec &vec ) const
{
  vec.set_size( solver->getSolution().n_rows );
  unsigned int extractCol = solver->getSolution().n_cols-1;
  for ( unsigned int i=0;i<vec.n_elem; i++ )
  {
    vec(i) = solver->getSolution()(i,extractCol).real();
  }
}

void WaveGuideFDSimulation::getExitField( arma::cx_vec &vec ) const
{
  vec.set_size( solver->getSolution().n_rows );
  unsigned int extractCol = solver->getSolution().n_cols-1;
  for ( unsigned int i=0;i<vec.n_elem; i++ )
  {
    vec(i) = solver->getSolution()(i,extractCol);
  }
}

void WaveGuideFDSimulation::computeFarField()
{
  // Extract the last column of the solution matrix
  arma::cx_vec exitField;
  getExitField( exitField );
  arma::cx_vec ft = arma::fft( exitField );
  farFieldModulus = new arma::vec( arma::abs(ft) );
}

void WaveGuideFDSimulation::saveFarField( const string &fname, unsigned int uid ) const
{
  arma::vec exitField;
  getExitField( exitField );
  hid_t file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  hsize_t dim = farFieldModulus->size();
  H5LTmake_dataset( file_id, "exitField", 1, &dim, H5T_NATIVE_DOUBLE, exitField.memptr());
  H5LTmake_dataset( file_id, "farField", 1, &dim, H5T_NATIVE_DOUBLE, farFieldModulus->memptr());
  H5LTset_attribute_double( file_id, "farField", "wavenumber", &wavenumber, 1);
  H5LTset_attribute_double( file_id, "farField", "gridspacing", &xDisc->step, 1);
  H5LTset_attribute_double( file_id, "exitField", "xmin", &xDisc->min, 1);
  H5LTset_attribute_double( file_id, "exitField", "xmax", &xDisc->max, 1);
  int intUID = uid;
  H5LTset_attribute_int( file_id, "farField", "uid", &intUID, 1);
  H5Fclose(file_id);
  clog << "Far field written to " << fname << endl;
}

cdouble WaveGuideFDSimulation::transverseBC( double z, Boundary_t bnd ) const
{
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

void WaveGuideFDSimulation::setBoundaryConditions( const ParaxialSource &source )
{
  src = &source;
  unsigned int Nx = nodeNumberTransverse();
  vector<cdouble> values(Nx, 1.0);
  for ( unsigned int i=0;i<Nx;i++ )
  {
    double x = xDisc->min + i*xDisc->step;
    values[i] = src->get( x, 0.0 );
  }
  solver->setLeftBC(&values[0]);
}
