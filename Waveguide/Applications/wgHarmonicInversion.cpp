#include <iostream>
#include "cladding.hpp"
#include "curvedWaveGuide2D.hpp"
#include "crankNicholson.hpp"
#include "controlFile.hpp"
#include "straightWG2D.hpp"
#include <complex>
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>
#include <harminv.h>
#include <jsoncpp/json/writer.h>
#include <H5Cpp.h>
#include <hdf5_hl.h>

using namespace std;
typedef complex<double> cdouble;

struct HarmInvRes
{
  vector<double> freq;
  vector<double> amplitude;
  vector<double> phase;
  vector<double> decay;
  double wcrd;
};

int main( int argc, char** argv )
{
  string ctlfname("data/singleCurvedWG702004.json");
  unsigned int nFreq = 100;
  double freqMin = 0.01;
  double freqMax = 10.0;

  try
  {
    ControlFile ctl;

    ctl.load(ctlfname);

    StraightWG2D wg;
    clog << "Loading solution... ";
    wg.init( ctl ); // Initialize the simulation
    clog << "done\n";

    double x0 = 0.0;
    double stepX = wg.transverseDiscretization().step;
    double width = wg.getWidth();
    double x = x0;

    vector<HarmInvRes> allRes;
    while ( x < x0+width )
    {
      clog << "Current x: " << x << " . Ends at: " << x0+width << endl;
      HarmInvRes results;
      vector<cdouble> res;
      results.wcrd = x;
      wg.extractField( x, res );

      // Put the imaginary part to zero (this is how harminv handles real signals)
      for ( unsigned int i=0;i<res.size();i++ )
      {
        res[i].imag(0.0);
      }

      harminv_data data = harminv_data_create( res.size(), &res[0], freqMin, freqMax, nFreq );

      harminv_solve( data );
      int foundFreq = harminv_get_num_freqs( data );

      if ( foundFreq == 0 )
      {
        cout << "Warning: Did not find any modes in the interval " << freqMin << "-" << freqMax << endl;
      }

      for ( unsigned int i=0;i<foundFreq;i++ )
      {
        results.freq.push_back(harminv_get_freq( data, i));
        results.decay.push_back(harminv_get_decay( data, i));
        cdouble amplitudeCompx;
        amplitudeCompx = harminv_get_amplitude(data, i);
        results.amplitude.push_back(abs(amplitudeCompx));
        results.phase.push_back( arg(amplitudeCompx) );
      }
      x += 20.0*stepX;
      allRes.push_back(results);
    }

    stringstream fname;
    fname << "data/harminvdata" << ctl.getUID() << ".h5";
    hid_t file_id = H5Fcreate(fname.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );;
    for ( unsigned int i=0;i<allRes.size();i++ )
    {
      stringstream freqname, decayname, ampname, phasename;
      freqname << "freq" << i;
      decayname << "decay" << i;
      ampname << "amplitude" << i;
      phasename << "phase" << i;
      hsize_t dim = allRes[i].freq.size();
      unsigned int rank = 1;

      H5LTmake_dataset( file_id, freqname.str().c_str(), rank, &dim, H5T_NATIVE_DOUBLE, &allRes[i].freq);
      H5LTmake_dataset( file_id, decayname.str().c_str(), rank, &dim, H5T_NATIVE_DOUBLE, &allRes[i].decay);
      H5LTmake_dataset( file_id, ampname.str().c_str(), rank, &dim, H5T_NATIVE_DOUBLE, &allRes[i].amplitude);
      H5LTmake_dataset( file_id, phasename.str().c_str(), rank, &dim, H5T_NATIVE_DOUBLE, &allRes[i].phase);

      // Set the attribute for the position
      H5LTset_attribute_double( file_id, freqname.str().c_str(), "x", &allRes[i].wcrd, 1);
      H5LTset_attribute_double( file_id, decayname.str().c_str(), "x", &allRes[i].wcrd, 1);
      H5LTset_attribute_double( file_id, ampname.str().c_str(), "x", &allRes[i].wcrd, 1);
      H5LTset_attribute_double( file_id, phasename.str().c_str(), "x", &allRes[i].wcrd, 1);
    }
    H5Fclose(file_id);
    clog << "Data written to " << fname.str() << endl;
  }
  catch ( exception &exc )
  {
    cout << exc.what() << endl;
    return 1;
  }
  return 0;
}
