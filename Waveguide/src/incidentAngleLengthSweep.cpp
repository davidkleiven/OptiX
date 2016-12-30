#include "incidentAngleLengthSweep.hpp"
#include <H5Cpp.h>
#include <hdf5_hl.h>

using namespace std;

void IncidentAngleLengthSweep::processStep( double z )
{
  if ( currentTheta >= nTheta )
  {
    cout << "Warning! Max number of angles reached!\n";
    return;
  }

  if ( ( z < nextZ ) || ( currentLengthIter >= Nlengths ) ) return;
  arma::vec farF;
  ff.result( wg.getSolver(), farF );

  // TODO: Move to set lmin
  //double z = wg.longitudinalDiscretization().max;
  //if ( z < Lmin )
  //{
    //throw (runtime_error("The minimum distance given is smaller than the maximum value in the discretization!"));
  //}

  /** Initialize the storage space */
  if ( data.size() == 0 )
  {
    data.resize( Nlengths );
    data[currentLengthIter].length = z;
  }

  /** Compute far field based on the at several positions */
  data[currentLengthIter].farField.set_size( farF.n_elem, nTheta );
  data[currentLengthIter].length = z;
  for ( unsigned int j=0;j<farF.n_elem;j++ )
  {
    data[currentLengthIter].farField(j, currentTheta) = farF(j);
  }

  double lengthStep = (zmax-Lmin)/static_cast<double>(Nlengths);
  nextZ += lengthStep;
  currentLengthIter++;
}

void IncidentAngleLengthSweep::proceedToNextAngle()
{
  nextZ = Lmin;
  currentLengthIter = 0;
  currentTheta++;
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
