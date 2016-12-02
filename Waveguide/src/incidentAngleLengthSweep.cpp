#include "incidentAngleLengthSweep.hpp"
#include <H5Cpp.h>
#include <hdf5_hl.h>

using namespace std;
void IncidentAngleLengthSweep::processFarField()
{
  if ( currentTheta >= nTheta )
  {
    cout << "Warning! Max number of angles reached!\n";
    return;
  }

  double z = wg.longitudinalDiscretization().max;
  if ( z < Lmin )
  {
    throw (runtime_error("The minimum distance given is smaller than the maximum value in the discretization!"));
  }
  double dz = (z-Lmin)/static_cast<double>(Nlengths);
  /** Initialize the storage space */
  if ( data.size() == 0 )
  {
    data.resize( Nlengths );
    for ( unsigned int i=0;i<Nlengths;i++ )
    {
      int n = angleIndx( 1.3*thetaMax );
      data[i].length = z-i*dz;
      data[i].farField.set_size( 2*n+1, nTheta );
    }
  }

  /** Compute far field based on the at several positions */
  for ( unsigned int i=0;i<Nlengths;i++ )
  {
    wg.computeFarField( fftSignalLength, z-i*dz );
    unsigned int n = data[i].farField.n_rows/2;
    for ( unsigned int j=0;j<2*n;j++ )
    {
      data[i].farField(j, currentTheta) = wg.getFarField()(fftSignalLength/2-n+j);
    }
  }
  currentTheta++; // Proceed to the next angle
}

void IncidentAngleLengthSweep::save( const string &fname ) const
{
  stringstream ss;
  if ( generateUID )
  {
    ss << fname << uid << ".h5";
  }
  else
  {
    ss << fname << ".h5";
  }

  hid_t file_id = H5Fcreate( ss.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hsize_t dim[2] = {data[0].farField.n_cols, data[0].farField.n_rows};
  double phiMin = angleDeg( -data[0].farField.n_rows/2);
  double phiMax = angleDeg( data[0].farField.n_rows/2);
  for ( unsigned int i=0;i<Nlengths;i++ )
  {
    stringstream dset;
    dset << "intensity" << i;
    H5LTmake_dataset( file_id, dset.str().c_str(), 2, dim, H5T_NATIVE_DOUBLE, data[i].farField.memptr());
    H5LTset_attribute_double( file_id, dset.str().c_str(), "thetaMin", &thetaMin, 1);
    H5LTset_attribute_double( file_id, dset.str().c_str(), "thetaMax", &thetaMax, 1);
    H5LTset_attribute_double( file_id, dset.str().c_str(), "phiMin", &phiMin, 1);
    H5LTset_attribute_double( file_id, dset.str().c_str(), "phiMax", &phiMax, 1);
    H5LTset_attribute_double( file_id, dset.str().c_str(), "position", &data[i].length, 1);
  }

  H5Fclose( file_id );
  clog << "Results written to " << ss.str() << endl;
}
