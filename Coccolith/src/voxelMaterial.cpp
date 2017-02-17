#include "voxelMaterial.hpp"
#include <iostream>
#include <stdexcept>

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

double VoxelMaterial::eps( const meep::vec &r )
{
  return 1.0;
}

double VoxelMaterial::conductivity( meep::component c, const meep::vec &r )
{
  return 0.0;
}
