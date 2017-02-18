#include "voxelMaterial.hpp"
#include <iostream>
#include <stdexcept>
#include <visa/visa.hpp>
#include <chrono>
#include <thread>

using namespace std;

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
  InfoFromFilename info;
  extractDimsFromFilename( fname, info );

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

  cout << "--------------------\n";
  cout << "Number of voxels outside material: " << numberOfZeros << endl;
  cout << "Number of voxels inside material: " << numberOfOnes << endl;
  cout << "--------------------\n";
}

void VoxelMaterial::slideThroughVoxelArray() const
{
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

  plots.get("Voxelvalues").clear();
  clog << "Showing YZ plane\n";
  for ( unsigned int i=0;i<voxels.n_rows;i++ )
  {
    arma::mat values;
    fillArmaMat( voxels.tube( i, 0, i, voxels.n_cols-1 ), values );
    plots.getActive().fillVertexArray( values );
    plots.show();
    this_thread::sleep_for( duration );
  }

  plots.get("Voxelvalues").clear();
  clog << "Showing XZ plane\n";
  for ( unsigned int i=0;i<voxels.n_cols;i++ )
  {
    arma::mat values;
    fillArmaMat( voxels.tube( 0, i, voxels.n_rows-1, i ), values );
    plots.getActive().fillVertexArray( values );
    plots.show();
    this_thread::sleep_for( duration );
  }
}

void VoxelMaterial::showProjections() const
{
  unsigned int holdForSec = 15;
  visa::WindowHandler plots;
  plots.addPlot( "XY-plane" );

  // Project in z-direction
  arma::mat values( voxels.n_rows, voxels.n_cols );
  values.fill(0.0);
  for ( unsigned int i=0;i<voxels.n_slices;i++ )
  {
    values +=  arma::conv_to<arma::mat>::from( voxels.slice(i) );
  }
  values /= voxels.n_slices;
  plots.getActive().fillVertexArray( values );

  plots.addPlot( "XZ-plane" );
  values.set_size( voxels.n_cols, voxels.n_slices );
  values.fill(0.0);
  for ( unsigned int i=0;i<voxels.n_rows;i++ )
  {
    arma::mat next;
    fillArmaMat( voxels.tube( i, 0, i, voxels.n_cols-1 ), next );
    values += next;
  }
  values /= voxels.n_rows;
  plots.get("XZ-plane").fillVertexArray( values );

  plots.addPlot( "YZ-plane" );
  values.set_size( voxels.n_rows, voxels.n_slices );
  values.fill(0.0);
  for ( unsigned int i=0;i<voxels.n_cols;i++ )
  {
    arma::mat next;
    fillArmaMat( voxels.tube( 0, i, voxels.n_rows-1, i ), next );
    values += next;
  }

  plots.get("YZ-plane").fillVertexArray( values );

  chrono::seconds duration(1);
  for ( unsigned int i=0;i<holdForSec;i++ )
  {
    plots.show();
    clog << "Closes in " << holdForSec-i << " seconds\r";
    this_thread::sleep_for( duration );
  }
  clog << endl;
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

////////////////////////////////////////////////////////////////////////////////
double CaCO3Cocco::eps( const meep::vec &r )
{
  if ( voxels( r.x(), r.y(), r.z() ) == 1 )
  {
    return epsilon;
  }
  return 1.0;
}

double CaCO3Cocco::conductivity( meep::component c, const meep::vec &r )
{
  return 0.0;
}
