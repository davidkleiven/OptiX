#include "paraxialSimulation.hpp"
#include "H5Cpp.h"
#include "solver2D.hpp"
#include "controlFile.hpp"
#include "crankNicholson.hpp"
#include "waveGuide.hpp"
#include "solver1D.hpp"
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include "paraxialSource.hpp"
#include "h5Attribute.hpp"
#include <limits>
#include <stdexcept>
#include <utility>

using namespace std;
const double PI = acos(-1.0);
ParaxialSimulation::ParaxialSimulation(const char* name):xDisc(new Disctretization), zDisc(new Disctretization), name(name){};

ParaxialSimulation::~ParaxialSimulation()
{
  if ( xDisc != NULL ) delete xDisc;
  if ( zDisc != NULL ) delete zDisc;
  if ( solverInitializedViaInit )
  {
    delete solver;
  }
  if ( farFieldModulus != NULL ) delete farFieldModulus;
  if ( file != NULL ) delete file;
}

void ParaxialSimulation::solve()
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

void ParaxialSimulation::setWaveLength( double lambda )
{
  wavenumber = 2.0*PI/lambda;
}

void ParaxialSimulation::setTransverseDiscretization( double xmin, double xmax, double step )
{
  xDisc->min = xmin;
  xDisc->max = xmax;
  xDisc->step = step;
}

void ParaxialSimulation::setLongitudinalDiscretization( double xmin, double xmax, double step )
{
  zDisc->min = xmin;
  zDisc->max = xmax;
  zDisc->step = step;
}

unsigned int ParaxialSimulation::nodeNumberTransverse() const
{
  return (xDisc->max - xDisc->min)/xDisc->step + 1.0;
}

unsigned int ParaxialSimulation::nodeNumberLongitudinal() const
{
  return ( zDisc->max - zDisc->min)/zDisc->step + 1.0;
}

void ParaxialSimulation::setSolver( Solver2D &solv )
{
  if ( solverInitializedViaInit )
  {
    throw ("You cannot set a new solver when it has been initialized via the init function");
  }

  solver = &solv;
  solver->setSimulator( *this );
}

void ParaxialSimulation::save( ControlFile &ctl )
{
  if ( solver == NULL )
  {
    throw ( runtime_error("No solver specified!\n"));
  }

  if ( src == NULL )
  {
    throw ( runtime_error("No source specified!\n") );
  }
  string fname = ctl.getFnameTemplate();
  string h5fname = fname+".h5";
  string jsonfname = fname+".json";
  vector<string> dsets;
  if ( file != NULL ) delete file;
  file = new H5::H5File( h5fname.c_str(), H5F_ACC_TRUNC );

  // Fill attrs
  commonAttributes.push_back( makeAttr("xmin",xDisc->min) );
  commonAttributes.push_back( makeAttr("xmax", xDisc->max) );
  commonAttributes.push_back( makeAttr("zmin", zDisc->min) );
  commonAttributes.push_back( makeAttr("zmax", zDisc->max) );
  int uid = ctl.getUID();
  commonAttributes.push_back( makeAttr("uid", uid) );

  //arma::abs(solver->getSolution()).save(h5fname.c_str(), arma::hdf5_binary);
  if ( saveColorPlot )
  {
    arma::mat absSol = arma::abs(solver->getSolution());
    hsize_t dim[2] = {absSol.n_rows, absSol.n_cols};
    int rank = 2;
    dsetnames.push_back("amplitude");
    saveArmaMat( absSol, dsets.back().c_str(), commonAttributes );

    solver->getPhase( absSol );
    dsetnames.push_back("phase");
    saveArmaMat( absSol, dsets.back().c_str(), commonAttributes );
  }

  saveFarField();

  Json::Value wginfo;
  Json::Value solverInfo;
  Json::Value sourceInfo;
  src->info( sourceInfo );
  ctl.get()["datafile"] = h5fname;
  ctl.get()["name"] = name;
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
}

void ParaxialSimulation::saveFarField()
{
  if ( farFieldModulus == NULL )
  {
    return;
  }

  arma::vec exitField;
  getExitField( exitField );
  arma::cx_vec exitFieldCmpl;
  getExitField( exitFieldCmpl );
  arma::vec exitAmpl = arma::abs( exitFieldCmpl );
  arma::vec exitPhase;
  exitPhase.set_size(exitAmpl.n_elem);
  for ( unsigned int i=0;i<exitAmpl.n_elem;i++ )
  {
    exitPhase(i) = arg( exitFieldCmpl(i) );
  }

  arma::vec zoomFar;
  extractFarField( zoomFar );
  hsize_t dim = zoomFar.n_elem;
  dsetnames.push_back("farField");
  saveArmaVec( zoomFar, dsetnames.back().c_str());
  H5::DataSet dset = file->openDataSet( dsetnames.back() );
  addAttribute( dset, "phiMin", farParam.phiMin );
  addAttribute( dset, "phiMax", farParam.phiMax );
  dim = exitField.n_elem;
  dsetnames.push_back("exitField");
  saveArmaVec( exitField, dsetnames.back().c_str() );
  dsetnames.push_back("exitIntensity");
  saveArmaVec( exitAmpl, dsetnames.back().c_str() );
  dsetnames.push_back("exitPhase");
  saveArmaVec( exitPhase, dsetnames.back().c_str() );
}

void ParaxialSimulation::getExitField( arma::vec &vec ) const
{
  vec.set_size( solver->getSolution().n_rows );
  unsigned int extractCol = solver->getSolution().n_cols-1;
  for ( unsigned int i=0;i<vec.n_elem; i++ )
  {
    vec(i) = solver->getSolution()(i,extractCol).real();
  }
}

void ParaxialSimulation::getExitField( arma::cx_vec &vec ) const
{
  vec.set_size( solver->getSolution().n_rows );
  unsigned int extractCol = solver->getSolution().n_cols-1;
  for ( unsigned int i=0;i<vec.n_elem; i++ )
  {
    vec(i) = solver->getSolution()(i,extractCol);
  }
}

void ParaxialSimulation::computeFarField( unsigned int signalLength )
{
  // Extract the last column, thus specify z to inf
  computeFarField( signalLength, numeric_limits<double>::max());
}

void ParaxialSimulation::computeFarField( unsigned int signalLength, double pos )
{

  // Extract the last column of the solution matrix
  arma::cx_vec exitField;
  arma::cx_vec paddedSignal;
  if ( pos > zDisc->max )
  {
    getExitField( exitField );
  }
  else
  {
    unsigned int ix, iz;
    closestIndex( 0.0, pos, ix, iz );
    iz = iz >= solver->getSolution().n_cols ? solver->getSolution().n_cols-1:iz;
    exitField.set_size( solver->getSolution().n_rows );
    for ( unsigned int ix=0;ix<solver->getSolution().n_rows;ix++ )
    {
      exitField(ix) = solver->getSolution()(ix,iz);
    }
  }
  if ( signalLength < exitField.n_elem )
  {
    paddedSignal = exitField;
  }
  else
  {
    paddedSignal.set_size(signalLength);
    paddedSignal.fill(farParam.padValue);
    // Fill in the signal on the center
    unsigned int start = signalLength/2 - exitField.n_elem/2;
    for ( unsigned int i=0;i<exitField.n_elem;i++ )
    {
      paddedSignal(start+i) = exitField(i);
    }
  }

  arma::cx_vec ft = arma::fft( paddedSignal );

  if ( farFieldModulus == NULL )
  {
    farFieldModulus = new arma::vec( arma::abs(ft)/sqrt(ft.n_elem) );
  }
  else
  {
    *farFieldModulus = arma::vec( arma::abs(ft)/sqrt(ft.n_elem) );
  }

  // Shift the FFT
  unsigned int N = farFieldModulus->n_elem;
  for ( unsigned int i=0;i<N/2;i++ )
  {
    double copy = (*farFieldModulus)(i);
    (*farFieldModulus)(i) = (*farFieldModulus)(i+N/2);
    (*farFieldModulus)(i+N/2) = copy;
  }
}

void ParaxialSimulation::computeFarField()
{
  computeFarField(0);
}

void ParaxialSimulation::closestIndex( double x, double z, unsigned int &ix, unsigned int &iz ) const
{
  ix = (x - xDisc->min)/xDisc->step;
  iz = (z - zDisc->min)/zDisc->step;
}

double ParaxialSimulation::getIntensity( double x, double z ) const
{
  // TODO: The interpolation in this routine does not work
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

  double intensity = pow( abs(solver->getSolution()(ix,iz)),2 )*(x2 - x)*(z2-z);
  intensity += pow( abs(solver->getSolution()(ix+1,iz)), 2 )*(x-x1)*(z2-z);
  intensity += pow( abs(solver->getSolution()(ix,iz+1)), 2 )*(z-z1)*(x2-x);
  intensity += pow( abs(solver->getSolution()(ix+1,iz+1)), 2 )*(z2-z)*(x-x1);
  return intensity/( (x2-x1)*(z2-z1) );
}

double ParaxialSimulation::getIntensity( unsigned int ix, unsigned int iz ) const
{
  return pow( abs( solver->getSolution()(ix,iz) ), 2 );
}

double ParaxialSimulation::getZ( unsigned int iz ) const
{
  return zDisc->min + iz*zDisc->step;
}

double ParaxialSimulation::getX( int ix ) const
{
  return xDisc->min + ix*xDisc->step;
}

void ParaxialSimulation::saveArmaMat( const arma::mat &matrix, const char* dsetname )
{
  saveArmaMat( matrix, dsetname, commonAttributes );
}

void ParaxialSimulation::saveArmaMat( const arma::mat &matrix, const char* dsetname, const vector<H5Attr> &attrs )
{
  // Create dataspace
  hsize_t fdim[2] = {matrix.n_rows, matrix.n_cols};
  H5::DataSpace dataspace( 2, fdim );

  // Create dataset
  H5::DataSet ds( file->createDataSet(dsetname, H5::PredType::NATIVE_DOUBLE, dataspace) );
  H5::DataSpace attribSpace(H5S_SCALAR);

  for ( unsigned int i=0;i<attrs.size();i++ )
  {
    H5::Attribute att = ds.createAttribute( attrs[i].name.c_str(), attrs[i].dtype, attribSpace );
    if ( attrs[i].dtype == H5::PredType::NATIVE_INT )
    {
      int value = attrs[i].value;
      att.write( H5::PredType::NATIVE_INT, &value );
    }
    double value = attrs[i].value;
    att.write( H5::PredType::NATIVE_DOUBLE, &value );
  }

  // Write to file
  ds.write( matrix.memptr(), H5::PredType::NATIVE_DOUBLE );
}

void ParaxialSimulation::saveArmaVec( const arma::vec &vec, const char* dsetname )
{
  saveArmaVec( vec, dsetname, commonAttributes );
}

void ParaxialSimulation::saveArmaVec( const arma::vec &vec, const char* dsetname, const vector<H5Attr> &attrs )
{
  // Create dataspace
  hsize_t fdim = vec.n_elem;
  H5::DataSpace dataspace( 1, &fdim );

  // Create dataset
  H5::DataSpace attribSpace(H5S_SCALAR);
  H5::DataSet ds( file->createDataSet(dsetname, H5::PredType::NATIVE_DOUBLE, dataspace) );

  for ( unsigned int i=0;i<attrs.size();i++ )
  {
    H5::Attribute att = ds.createAttribute( attrs[i].name.c_str(), attrs[i].dtype, attribSpace );
    if ( attrs[i].dtype == H5::PredType::NATIVE_INT )
    {
      int value = attrs[i].value;
      att.write( H5::PredType::NATIVE_INT, &value );
    }
    double value = attrs[i].value;
    att.write( H5::PredType::NATIVE_DOUBLE, &value );
  }

  // Write to file
  ds.write( vec.memptr(), H5::PredType::NATIVE_DOUBLE );
}

void ParaxialSimulation::saveVec( const vector<double> &vec, const char* dsetname )
{
  // Create dataspace
  hsize_t fdim = vec.size();
  H5::DataSpace dataspace( 1, &fdim );

  // Create dataset
  H5::DataSet dataset( file->createDataSet(dsetname, H5::PredType::NATIVE_DOUBLE, dataspace) );

  // Write to file
  dataset.write( &vec[0], H5::PredType::NATIVE_DOUBLE );
}

void ParaxialSimulation::addAttribute( H5::DataSet &ds, const char* name, double value )
{
  H5::DataSpace attribSpace(H5S_SCALAR);
  H5::Attribute att = ds.createAttribute( name, H5::PredType::NATIVE_DOUBLE, attribSpace );
  att.write( H5::PredType::NATIVE_DOUBLE, &value );
}

void ParaxialSimulation::addAttribute( H5::DataSet &ds, const char* name, int value )
{
  H5::DataSpace attribSpace(H5S_SCALAR);
  H5::Attribute att = ds.createAttribute( name, H5::PredType::NATIVE_INT, attribSpace );
  att.write( H5::PredType::NATIVE_INT, &value );
}

void ParaxialSimulation::setBoundaryConditions( const ParaxialSource &source )
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

double ParaxialSimulation::getEnergy() const
{
  return 12.398*100.0/getWavelength();
}

double ParaxialSimulation::getWavelength() const
{
  return 2.0*PI/wavenumber;
}

void ParaxialSimulation::extractFarField( arma::vec &newFarField ) const
{
  if ( farFieldModulus == NULL )
  {
    throw (runtime_error("ExtractFarField: Far field is not computed!"));
  }
  unsigned int nmin = farFieldAngleToIndx( farParam.phiMin );
  unsigned int nmax = farFieldAngleToIndx( farParam.phiMax );
  newFarField = farFieldModulus->subvec(nmin, nmax);
}

unsigned int ParaxialSimulation::farFieldAngleToIndx( double angle ) const
{
  if ( farFieldModulus == NULL )
  {
    throw (runtime_error("farFieldAngleToIndx: Far field is not computed!"));
  }
  double angleRad = angle*PI/180.0;
  int n = farFieldModulus->n_elem*wavenumber*xDisc->step*sin(angleRad)/(2.0*PI);
  if ( abs(n) >= farFieldModulus->n_elem/2 )
  {
    if ( angle > 0.0 ) return farFieldModulus->n_elem-1;
    return 0;
  }
  return n+static_cast<int>(farFieldModulus->n_elem)/2;
}

void ParaxialSimulation::setFarFieldAngleRange( double phiMin, double phiMax )
{
  farParam.phiMin = phiMin;
  farParam.phiMax = phiMax;
}
