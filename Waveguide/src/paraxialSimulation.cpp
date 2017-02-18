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
#include "arraySource.hpp"
#include "h5Attribute.hpp"
#include "hdf5DataspaceCreator.hpp"
#include <limits>
#include <stdexcept>
#include <utility>
#include <sstream>

using namespace std;
const double PI = acos(-1.0);
ParaxialSimulation::ParaxialSimulation(const char* name):xDisc(new Disctretization), zDisc(new Disctretization),
yDisc(new Disctretization), name(name){};

ParaxialSimulation::~ParaxialSimulation()
{
  if ( xDisc != NULL ) delete xDisc;
  if ( zDisc != NULL ) delete zDisc;
  if ( yDisc != NULL ) delete yDisc;

  if ( solverInitializedViaInit )
  {
    delete solver;
  }
  if ( farFieldModulus != NULL ) delete farFieldModulus;
  if ( file != NULL ) delete file;
  if ( maingroup != NULL ) delete maingroup;
}

void ParaxialSimulation::solve()
{
  solver->solve();
}

void ParaxialSimulation::step()
{
  solver->step();
}

void ParaxialSimulation::reset()
{
  solver->reset();
}

void ParaxialSimulation::verifySolverReady() const
{
  if ( solver == NULL )
  {
    throw ( runtime_error("No solver specified!"));
  }

  if ( src == NULL )
  {
    throw ( runtime_error("No source specified!"));
  }
}

void ParaxialSimulation::setWaveLength( double lambda )
{
  wavenumber = 2.0*PI/lambda;
}

void ParaxialSimulation::setTransverseDiscretization( double xmin, double xmax, double step )
{
  setTransverseDiscretization( xmin, xmax, step, 1 );
}

void ParaxialSimulation::setTransverseDiscretization( double xmin, double xmax, double step, unsigned int downsamplingRatio )
{
  xDisc->min = xmin;
  xDisc->max = xmax;
  xDisc->step = step;
  xDisc->downsamplingRatio = downsamplingRatio;
  yDisc->downsamplingRatio = xDisc->downsamplingRatio;
}

void ParaxialSimulation::setVerticalDiscretization( double ymin, double ymax, double step )
{
  yDisc->min = ymin;
  yDisc->max = ymax;
  yDisc->step = step;
}

void ParaxialSimulation::setLongitudinalDiscretization( double zmin, double zmax, double step )
{
  setLongitudinalDiscretization( zmin, zmax, step, 1 );
}

void ParaxialSimulation::setLongitudinalDiscretization( double zmin, double zmax, double step, unsigned int downsamplingRatio )
{
  zDisc->min = zmin;
  zDisc->max = zmax;
  zDisc->step = step;
  zDisc->downsamplingRatio = downsamplingRatio;
}

unsigned int ParaxialSimulation::nodeNumberTransverse() const
{
  return (xDisc->max - xDisc->min)/xDisc->step + 1.0;
}

unsigned int ParaxialSimulation::nodeNumberVertical() const
{
  return (yDisc->max - yDisc->min)/yDisc->step + 1.0;
}

unsigned int ParaxialSimulation::nodeNumberLongitudinal() const
{
  return ( zDisc->max - zDisc->min)/zDisc->step + 1.0;
}

void ParaxialSimulation::setSolver( Solver &solv )
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
  maingroup = new H5::Group( file->createGroup(groupname+"/") );
  uid = ctl.getUID();
  setGroupAttributes();

  // Save all results from all the post processing modules
  for ( unsigned int i=0;i<postProcess.size();i++ )
  {
    vector<H5Attr> attrib;
    arma::vec res1D;
    arma::mat res2D;
    arma::cube res3D;

    // Get the result. No two of these functions will always be empty
    postProcess[i]->addAttrib( attrib );
    switch ( postProcess[i]->getReturnType( *solver ) )
    {
      case ( post::PostProcessingModule::ReturnType_t::vector1D ):
        postProcess[i]->result( *solver, res1D );
        saveArray( res1D, postProcess[i]->getName().c_str(), attrib );
        break;
      case ( post::PostProcessingModule::ReturnType_t::matrix2D ):
        postProcess[i]->result( *solver, res2D );
        saveArray( res2D, postProcess[i]->getName().c_str(), attrib );
        break;
      case ( post::PostProcessingModule::ReturnType_t::cube3D ):
        postProcess[i]->result( *solver, res3D );
        saveArray( res3D, postProcess[i]->getName().c_str(), attrib );
        break;
    }
    clog << "Dataset " << postProcess[i]->getName() << " added to HDF5 file\n";
  }

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

double ParaxialSimulation::getY( int iy ) const
{
  return yDisc->min + iy*yDisc->step;
}

template <class arrayType>
void ParaxialSimulation::saveArray( arrayType &matrix, const char* dsetname )
{
  vector<H5Attr> dummy;
  saveArray( matrix, dsetname, dummy );
}

template <class arrayType>
void ParaxialSimulation::saveArray( arrayType &matrix, const char* dsetname, const vector<H5Attr> &attrs )
{
  // Create dataspace
  hsize_t fdim[3];
  DataspaceCreator<arrayType> dsinfo;
  dsinfo.setDims( matrix, fdim );

  H5::DataSpace dataspace( dsinfo.rank, fdim );

  string name(groupname);
  name += dsetname;
  // Create dataset
  H5::DataSet ds( file->createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace) );
  H5::DataSpace attribSpace(H5S_SCALAR);
  for ( unsigned int i=0;i<attrs.size();i++ )
  {
    H5::Attribute att = ds.createAttribute( attrs[i].name.c_str(), attrs[i].dtype, attribSpace );
    if ( attrs[i].dtype == H5::PredType::NATIVE_INT )
    {
      int value = attrs[i].value;
      att.write( H5::PredType::NATIVE_INT, &value );
    }
    else
    {
      double value = attrs[i].value;
      att.write( H5::PredType::NATIVE_DOUBLE, &value );
    }
  }

  // Write to file
  ds.write( matrix.memptr(), H5::PredType::NATIVE_DOUBLE );
}

// Pre-fine allowed template types
template void ParaxialSimulation::saveArray<arma::vec>( arma::vec &matrix, const char* dsetname, const vector<H5Attr> &attrs );
template void ParaxialSimulation::saveArray<arma::mat>( arma::mat &matrix, const char* dsetname, const vector<H5Attr> &attrs );
template void ParaxialSimulation::saveArray<arma::cube>( arma::cube &matrix, const char* dsetname, const vector<H5Attr> &attrs );
template void ParaxialSimulation::saveArray<arma::vec>( arma::vec &matrix, const char* dsetname );
template void ParaxialSimulation::saveArray<arma::mat>( arma::mat &matrix, const char* dsetname );
template void ParaxialSimulation::saveArray<arma::cube>( arma::cube &m, const char* dsetname );

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

  typedef ParaxialSource::Dim_t Dim_t;
  switch ( source.getDim() )
  {
    case Dim_t::TWO_D:
    {
      arma::cx_vec values(Nx, 1.0);
      for ( unsigned int i=0;i<Nx;i++ )
      {
        double x = getX(i);
        values[i] = src->get( x, 0.0 );
      }
      solver->setInitialConditions( values );
      break;
    }
    case Dim_t::THREE_D:
    {
      arma::cx_mat values(Nx, Nx);
      for ( unsigned int i=0;i<Nx;i++ )
      {
        double x = getX(i);
        for ( unsigned int j=0;j<Nx;j++ )
        {
          double y = getY(j);
          values(j,i) = src->get(x,y,0.0);
        }
      }
      solver->setInitialConditions( values );
    }
  }
}

void ParaxialSimulation::setBoundaryConditions( const ArraySource &source )
{
  src = &source;
  if ( source.getVec().n_elem != nodeNumberTransverse() )
  {
    stringstream msg;
    msg << "The size of the given array does not match the number of transverse nodes! ";
    msg << "Required size: " << nodeNumberTransverse();
    msg << " Given size: " << source.getVec().n_elem;
    throw ( runtime_error( msg.str() ) );
  }
  solver->setInitialConditions( source.getVec() );
}

double ParaxialSimulation::getEnergy() const
{
  return 12.398*100.0/getWavelength();
}

double ParaxialSimulation::getWavelength() const
{
  return 2.0*PI/wavenumber;
}

void ParaxialSimulation::setFarFieldAngleRange( double phiMin, double phiMax )
{
  farParam.phiMin = phiMin;
  farParam.phiMax = phiMax;
}

ParaxialSimulation& ParaxialSimulation::operator << ( post::PostProcessingModule &module )
{
  postProcess.push_back( &module );
  return *this;
}

ParaxialSimulation& ParaxialSimulation::operator << (post::FarField &ff )
{
  ff.linkParaxialSim( *this );
  postProcess.push_back( &ff );
  return *this;
}

void ParaxialSimulation::setGroupAttributes()
{
  if ( maingroup == NULL ) return;

  // Create a dataspace
  H5::DataSpace attribSpace(H5S_SCALAR);

  // Write first attribute
  H5::Attribute att = maingroup->createAttribute( "xmin", H5::PredType::NATIVE_DOUBLE, attribSpace );
  att.write( H5::PredType::NATIVE_DOUBLE, &xDisc->min );

  att = maingroup->createAttribute( "xmax", H5::PredType::NATIVE_DOUBLE, attribSpace );
  att.write( H5::PredType::NATIVE_DOUBLE, &xDisc->max );

  att = maingroup->createAttribute( "dx", H5::PredType::NATIVE_DOUBLE, attribSpace );
  att.write( H5::PredType::NATIVE_DOUBLE, &xDisc->step );

  att = maingroup->createAttribute( "zmin", H5::PredType::NATIVE_DOUBLE, attribSpace );
  att.write( H5::PredType::NATIVE_DOUBLE, &zDisc->min );

  att = maingroup->createAttribute( "zmax", H5::PredType::NATIVE_DOUBLE, attribSpace );
  att.write( H5::PredType::NATIVE_DOUBLE, &zDisc->max );

  att = maingroup->createAttribute( "dz", H5::PredType::NATIVE_DOUBLE, attribSpace );
  att.write( H5::PredType::NATIVE_DOUBLE, &zDisc->step );

  att = maingroup->createAttribute( "ymin", H5::PredType::NATIVE_DOUBLE, attribSpace );
  att.write( H5::PredType::NATIVE_DOUBLE, &yDisc->min );

  att = maingroup->createAttribute( "ymax", H5::PredType::NATIVE_DOUBLE, attribSpace );
  att.write( H5::PredType::NATIVE_DOUBLE, &yDisc->max );

  att = maingroup->createAttribute( "dy", H5::PredType::NATIVE_DOUBLE, attribSpace );
  att.write( H5::PredType::NATIVE_DOUBLE, &yDisc->step );

  att = maingroup->createAttribute( "uid", H5::PredType::NATIVE_INT, attribSpace );
  att.write( H5::PredType::NATIVE_INT, &uid );

  att = maingroup->createAttribute( "wavenumber", H5::PredType::NATIVE_DOUBLE, attribSpace );
  att.write( H5::PredType::NATIVE_DOUBLE, &wavenumber );
}
