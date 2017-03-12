#include "voxelMaterial.hpp"
#include "config.h"
#include <iostream>
#include <stdexcept>
#ifdef HAVE_LIB_VISA
  #include <visa/visa.hpp>
#endif
#include <chrono>
#include <thread>

using namespace std;

arma::Cube<unsigned char> VoxelMaterial::voxels;
bool VoxelMaterial::materialIsLoaded = false;
bool VoxelMaterial::referenceRun = false;
double VoxelMaterial::vxsize = 1.0;
DomainSize VoxelMaterial::domain;

void VoxelMaterial::extractDimsFromFilename( const string &fname, InfoFromFilename &info )
{
  unsigned int values[4];
  unsigned int pos = 0;
  for ( unsigned int i=0;i<4;i++ )
  {
    pos = fname.find("_", pos+1);

    if ( pos == string::npos )
    {
      throw ( runtime_error("Filename does not contain 4 numbers separated with underscore!") );
    }

    stringstream ss;
    ss << fname.substr(pos+1,3);
    ss >> values[i];
  }
  info.voxelsize = static_cast<double>(values[0])/10.0;
  info.Nx = values[1];
  info.Ny = values[2];
  info.Nz = values[3];
}

void VoxelMaterial::loadRaw( const char* fname )
{
  string strfname(fname);
  loadRaw(strfname);
}

void VoxelMaterial::loadRaw( const string &fname )
{
  if ( materialIsLoaded )
  {
    clog << "Warning! The voxels have already been loaded. This will change the behaviour of all voxel materials\n";
  }
  InfoFromFilename info;
  extractDimsFromFilename( fname, info );
  vxsize = info.voxelsize;

  voxels.load( fname, arma::raw_binary );

  // Check that the total size adds up
  unsigned int expectedSize = info.Nx*info.Ny*info.Nz;
  if ( expectedSize != voxels.n_elem )
  {
    stringstream ss;
    ss << "Read wrong size from file. Expected: " << expectedSize << ". Read: " << voxels.n_elem << ".";
    throw( runtime_error(ss.str()) );
  }

  voxels.reshape( info.Nx, info.Ny, info.Nz );
  applyThreshold();
  showStatistics();
  materialIsLoaded = true;
}

void VoxelMaterial::applyThreshold()
{
  unsigned char maxval = voxels.max();
  for ( unsigned int k=0;k<voxels.n_slices;k++ )
    for ( unsigned int i=0;i<voxels.n_cols;i++ )
      for ( unsigned int j=0;j<voxels.n_rows;j++ )
      {
        if ( voxels(j,i,k) > 0.5*maxval )
        {
          voxels(j,i,k) = 1;
        }
        else
        {
          voxels(j,i,k) = 0;
        }
      }
}

void VoxelMaterial::showStatistics() const
{
  unsigned int numberOfOnes = 0;
  unsigned int numberOfZeros = 0;

  for ( unsigned int k=0;k<voxels.n_slices;k++ )
    for ( unsigned int i=0;i<voxels.n_cols;i++ )
      for ( unsigned int j=0;j<voxels.n_rows;j++ )
      {
        if ( voxels(j,i,k) == 0 )
        {
          numberOfZeros++;
        }
        else
        {
          numberOfOnes++;
        }
      }

  if ( meep::am_master() )
  {
    cout << "---------------------------------------------------------------\n";
    cout << "Number of voxels outside material: " << numberOfZeros << endl;
    cout << "Number of voxels inside material: " << numberOfOnes << endl;
    cout << "---------------------------------------------------------------\n";
  }
}

void VoxelMaterial::slideThroughVoxelArray() const
{
  #ifdef HAVE_LIB_VISA
    clog << "Showing XY plane\n";
    visa::WindowHandler plots;
    plots.addPlot( "Voxelvalues" );
    chrono::milliseconds duration(50);
    for ( unsigned int i=0;i<voxels.n_slices;i++ )
    {
      arma::mat values = arma::conv_to<arma::mat>::from( voxels.slice(i) );
      plots.getActive().fillVertexArray( values );
      plots.show();
      this_thread::sleep_for( duration );
    }

    plots.clear();
    clog << "Showing YZ plane\n";
    for ( unsigned int i=0;i<voxels.n_rows;i++ )
    {
      arma::mat values;
      fillArmaMat( voxels.tube( i, 0, i, voxels.n_cols-1 ), values );
      plots.getActive().fillVertexArray( values );
      plots.show();
      this_thread::sleep_for( duration );
    }

    clog << "Showing XZ plane\n";
    for ( unsigned int i=0;i<voxels.n_cols;i++ )
    {
      arma::mat values;
      fillArmaMat( voxels.tube( 0, i, voxels.n_rows-1, i ), values );
      plots.getActive().fillVertexArray( values );
      plots.show();
      this_thread::sleep_for( duration );
    }
  #else
    clog << "Compiled without libVISA!\n";
  #endif
}

void VoxelMaterial::showProjections() const
{
  #ifdef HAVE_LIB_VISA
    unsigned int holdForSec = 15;
    visa::WindowHandler plots;
    plots.setLayout(2,2);
    plots.addPlot( "XY-plane" );
    plots.get( "XY-plane" ).setCmap( visa::Colormaps::Colormap_t::GREYSCALE );

    // Project in z-direction
    arma::mat values( voxels.n_rows, voxels.n_cols );
    projectionXY( values );
    plots.getActive().setImg( values );

    plots.addPlot( "XZ-plane" );
    plots.get( "XZ-plane" ).setCmap( visa::Colormaps::Colormap_t::GREYSCALE );
    values.set_size( voxels.n_cols, voxels.n_slices );
    projectionXZ( values );
    plots.get("XZ-plane").setImg( values );

    plots.addPlot( "YZ-plane" );
    plots.get( "YZ-plane" ).setCmap( visa::Colormaps::Colormap_t::GREYSCALE );
    values.set_size( voxels.n_rows, voxels.n_slices );

    projectionYZ( values );
    plots.get("YZ-plane").setImg( values );

    chrono::seconds duration(1);
    for ( unsigned int i=0;i<holdForSec;i++ )
    {
      plots.show();
      clog << "Closes in " << holdForSec-i << " seconds\r";
      this_thread::sleep_for( duration );
    }
    clog << endl;
  #else
    clog << "Compiled without libvisa!\n";
  #endif
}

void VoxelMaterial::fillArmaMat( const arma::Mat<unsigned char> &values, arma::mat &matrix )
{
  matrix.set_size( values.n_rows, values.n_cols );
  for ( unsigned int i=0;i<values.n_cols;i++ )
  {
    for ( unsigned int j=0;j<values.n_rows;j++ )
    {
      matrix(j,i) = values(j,i);
    }
  }
}

void VoxelMaterial::projectionXY( arma::mat &matrix ) const
{
  matrix.set_size( voxels.n_rows, voxels.n_cols );
  matrix.fill(0.0);
  for ( unsigned int i=0;i<voxels.n_slices;i++ )
  {
    matrix +=  arma::conv_to<arma::mat>::from( voxels.slice(i) );
  }
  matrix /= voxels.n_slices;
}

void VoxelMaterial::projectionXZ( arma::mat &matrix ) const
{
  matrix.fill(0.0);
  for ( unsigned int i=0;i<voxels.n_rows;i++ )
  {
    arma::mat next;
    fillArmaMat( voxels.tube( i, 0, i, voxels.n_cols-1 ), next );
    matrix += next;
  }
  matrix /= voxels.n_rows;
}

void VoxelMaterial::projectionYZ( arma::mat &matrix ) const
{
  matrix.fill(0.0);
  for ( unsigned int i=0;i<voxels.n_cols;i++ )
  {
    arma::mat next;
    fillArmaMat( voxels.tube( 0, i, voxels.n_rows-1, i ), next );
    matrix += next;
  }
  matrix /= voxels.n_cols;
}

void VoxelMaterial::setDomainSize( const meep::grid_volume &gvol, double PMLThick )
{
  // The bounding box of this domain corresponds to the data read from the file
  // There is one PML layer on all sides, this this is added/subtracted from the
  // max/min values in the simulation volume
  domain.xmin = gvol.xmin()+PMLThick;
  domain.xmax = gvol.xmax()-PMLThick;
  domain.ymin = gvol.ymin()+PMLThick;
  domain.ymax = gvol.ymax()-PMLThick;
  domain.zmin = gvol.zmin()+PMLThick;
  domain.zmax = gvol.zmax()-PMLThick;
}

void VoxelMaterial::meepVecToIndx( const meep::vec &r, unsigned int indx[3] )
{
    assert( isInsideDomain(r) );
    indx[0] = ( r.x() - domain.xmin )*(voxels.n_rows-1)/( domain.xmax - domain.xmin );
    indx[1] = (r.y() - domain.ymin)*(voxels.n_cols-1)/( domain.ymax-domain.ymin );
    indx[2] = ( r.z() - domain.zmin )*( voxels.n_slices-1)/( domain.zmax - domain.zmin );
}

bool VoxelMaterial::isInsideDomain( const meep::vec &r )
{
  return ( r.x() > domain.xmin ) && ( r.x() < domain.xmax ) && \
         ( r.y() > domain.ymin ) && ( r.y() < domain.ymax ) && \
         ( r.z() > domain.zmin ) && ( r.z() < domain.zmax );
}

void VoxelMaterial::boundingBox( unsigned int lowerCrn[3], unsigned int upperCrn[3] ) const
{
  lowerCrn[0] = voxels.n_rows-1;
  lowerCrn[1] = voxels.n_cols-1;
  lowerCrn[2] = voxels.n_slices-1;
  upperCrn[0] = 0;
  upperCrn[1] = 0;
  upperCrn[2] = 0;
  for ( unsigned int z=0;z<voxels.n_slices;z++ )
  for ( unsigned int y=0;y<voxels.n_cols;y++ )
  for ( unsigned int x=0;x<voxels.n_rows;x++ )
  {
    if ( voxels(x,y,z) == 1 )
    {
      if ( x < lowerCrn[0] ) lowerCrn[0] = x;
      if ( x > upperCrn[0] ) upperCrn[0] = x;
      if ( y < lowerCrn[1] ) lowerCrn[1] = y;
      if ( y > upperCrn[1] ) upperCrn[1] = y;
      if ( z < lowerCrn[2] ) lowerCrn[2] = z;
      if ( z > upperCrn[2] ) upperCrn[2] = z;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
double CaCO3Cocco::eps( const meep::vec &r )
{
  if ( referenceRun )
  {
    return 1.0;
  }
  else if ( !isInsideDomain(r) )
  {
    return 1.0;
  }

  unsigned int indx[3];
  meepVecToIndx( r, indx );

  if ( voxels( indx[0], indx[1], indx[2] ) == 1 )
  {
    return epsilon;
  }
  return 1.0;
}

double CaCO3Cocco::conductivity( meep::component c, const meep::vec &r )
{
  return 0.0;
}

double CaCO3Cocco::chi1p1( meep::field_type ft, const meep::vec &r )
{
  if ( ft != meep::E_stuff ) return 1.0;

  return eps(r);
}

////////////////////////////////////////////////////////////////////////////////
void VoxelSusceptibility::sigma_row( meep::component c, double sigrow[3], const meep::vec &r )
{
  sigrow[0] = sigrow[1] = sigrow[2] = 0.0;
  sigrow[meep::component_index(c)] = f(r);
}

double VoxelSusceptibility::f( const meep::vec &r )
{
  if ( isReferenceRun() || !isInsideDomain(r)) return referenceReturnVal;
  unsigned int indx[3];
  meepVecToIndx( r, indx );

  if ( voxels( indx[0], indx[1], indx[2] ) == 1 )
  {
    return parameter;
  }
  return referenceReturnVal;
}

////////////////////////////////////////////////////////////////////////////////
