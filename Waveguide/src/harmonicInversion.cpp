#include "harmonicInversion.hpp"
#include <harminv.h>
#include <utility>
#include <iostream>
#include <H5Cpp.h>
#include <hdf5_hl.h>

using namespace std;

void HarmonicInversion::solve( vector<cdouble> &values )
{
  harminv_data data = harminv_data_create( values.size(), &values[0], freqMin, freqMax, nFreq );
  harminv_solve( data );
  int foundFreq = harminv_get_num_freqs( data );

  if ( foundFreq == 0 )
  {
    cout << "Warning: Did not find any modes in the interval " << freqMin << "-" << freqMax << endl;
  }

  HarmInvRes res;
  for ( unsigned int i=0;i<values.size();i++ )
  {
    res.data.push_back( values[i].real() );
  }

  for ( unsigned int i=0;i<foundFreq;i++ )
  {
    res.freq.push_back(harminv_get_freq( data, i));
    res.decay.push_back(harminv_get_decay( data, i));
    cdouble amplitudeCompx;
    amplitudeCompx = harminv_get_amplitude(data, i);
    res.amplitude.push_back(abs(amplitudeCompx));
    res.phase.push_back( arg(amplitudeCompx) );
  }
  results.push_back(res);
}

void HarmonicInversion::setFreq( double frMin, double frMax, unsigned int nfr )
{
  freqMin = frMin;
  freqMax = frMax;
  nFreq = nfr;
}

void HarmonicInversion::addAttribute( const char* name, double value, unsigned int dset )
{
  if ( dset >= results.size() )
  {
    cout << "Dataset does not exists\n";
    return;
  }
  results[dset].attrib.insert( pair<string,double>(name,value) );
}

void HarmonicInversion::save( const char* fname ) const
{
  hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  for ( unsigned int i=0;i<results.size();i++ )
  {
    stringstream freqname, decayname, ampname, phasename, dataname;
    freqname << "freq" << i;
    decayname << "decay" << i;
    ampname << "amplitude" << i;
    phasename << "phase" << i;
    dataname << "transformedData" << i;
    hsize_t dim = results[i].freq.size();
    unsigned int rank = 1;

    H5LTmake_dataset( file_id, freqname.str().c_str(), rank, &dim, H5T_NATIVE_DOUBLE, &results[i].freq[0]);
    H5LTmake_dataset( file_id, decayname.str().c_str(), rank, &dim, H5T_NATIVE_DOUBLE, &results[i].decay[0]);
    H5LTmake_dataset( file_id, ampname.str().c_str(), rank, &dim, H5T_NATIVE_DOUBLE, &results[i].amplitude[0]);
    H5LTmake_dataset( file_id, phasename.str().c_str(), rank, &dim, H5T_NATIVE_DOUBLE, &results[i].phase[0]);
    dim = results[i].data.size();
    H5LTmake_dataset( file_id, dataname.str().c_str(), rank, &dim, H5T_NATIVE_DOUBLE, &results[i].data[0]);

    for ( auto iter=results[i].attrib.begin(); iter != results[i].attrib.end(); ++iter )
    {
      // Set the attribute for the position
      H5LTset_attribute_double( file_id, freqname.str().c_str(), iter->first.c_str(), &iter->second, 1);
      H5LTset_attribute_double( file_id, decayname.str().c_str(), iter->first.c_str(), &iter->second, 1);
      H5LTset_attribute_double( file_id, ampname.str().c_str(), iter->first.c_str(), &iter->second, 1);
      H5LTset_attribute_double( file_id, phasename.str().c_str(), iter->first.c_str(), &iter->second, 1);
      H5LTset_attribute_double( file_id, dataname.str().c_str(), iter->first.c_str(), &iter->second, 1);
    }
  }
  H5Fclose(file_id);
  clog << "Data written to " << fname << endl;
}
