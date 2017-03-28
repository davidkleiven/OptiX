#include "voxelMaterial2D.hpp"
#include <stdexcept>
#include <cassert>

using namespace std;

arma::Mat<unsigned char> VoxelMaterial2D::voxels;

void VoxelMaterial2D::extractInfoFromFilename( const string &fname )
{
  unsigned int values[3];
  unsigned int pos = 0;
  for ( unsigned int i=0;i<3;i++ )
  {
    pos = fname.find("_", pos+1);
    if ( pos == string::npos )
    {
      throw (runtime_error("Wrong filename format. Expected <somefilename>_resolution_Nx_Ny.raw"));
    }
    stringstream ss;
    unsigned int length = 0;

    if ( !isdigit(fname[pos+1]) )
    {
      throw( runtime_error("The first number after the underscore is not a digit!") );
    }
    for (unsigned int j=pos+1;j<fname.length();j++ )
    {
      if ( !isdigit(fname[j]) ) break;
      length++;
    }

    ss << fname.substr( pos+1, length );
    ss >> values[i];
  }

  info.voxelsize = values[0]/10.0;
  info.Nx = values[1];
  info.Ny = values[2];
}

void VoxelMaterial2D::loadRaw( const string &fname )
{
  voxels.load( fname.c_str(), arma::raw_binary );
  extractInfoFromFilename( fname );

  if ( info.Nx*info.Ny != voxels.n_elem )
  {
    throw( runtime_error("The dimensions from the filename, does not match the number of elements read!"));
  }
  voxels.reshape( info.Nx, info.Ny );
  applyThreshold();
}

void VoxelMaterial2D::loadRaw( const char* fname )
{
  string fnameCppString(fname);
  loadRaw(fnameCppString);
}

void VoxelMaterial2D::applyThreshold()
{
  for ( unsigned int i=0;i<voxels.n_cols;i++ )
  {
    for ( unsigned int j=0;j<voxels.n_rows;j++ )
    {
      voxels(j,i) = voxels(j,i) > threshold ? 1:0;
    }
  }
}

bool VoxelMaterial2D::isOutsideComputationalDomain( const meep::vec &r ) const
{
  return (r.x() > domain.xmax ) || (r.x() < domain.xmin ) || ( r.y() < domain.ymin ) || ( r.y() > domain.ymax );
}

unsigned char VoxelMaterial2D::get( const meep::vec &r ) const
{
  assert( !isOutsideComputationalDomain(r) );
  unsigned int x = Nx()*( r.x()-domain.xmin)/(domain.xmax-domain.xmin);
  unsigned int y = Ny()*( r.y()-domain.ymin)/(domain.ymax-domain.ymin);
  return get(x,y);
}

bool VoxelMaterial2D::isInside( const meep::vec& r ) const
{
  if ( isOutsideComputationalDomain(r) ) return false;
  return get(r) == 1;
}

// ===================== DISPERSIVE 2D VOXEL MATERIAL ==========================
double Voxel2DSusceptibility::f( const meep::vec &r )
{
  if ( isInside(r) ) return sigma;
  return defaultParam;
}

void Voxel2DSusceptibility::sigma_row( meep::component c, double sigrow[3], const meep::vec &r )
{
  sigrow[0] = sigrow[1] = sigrow[2] = 0.0;
  sigrow[meep::component_index(c)] = f(r);
}
