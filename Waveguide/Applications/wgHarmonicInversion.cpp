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
#include <armadillo>

using namespace std;
typedef complex<double> cdouble;

struct HarmInvRes
{
  vector<double> freq;
  vector<double> amplitude;
  vector<double> phase;
  vector<double> decay;
  vector<double> data;
  double wcrd;
};

int main( int argc, char** argv )
{
  string ctlfname("data/singleCurvedWG966449.json");
  unsigned int nFreq = 10;
  double freqMin = 1E-6;
  double freqMax = 0.005; // Nyquist
  bool useFFT = false;

  try
  {
    ControlFile ctl;

    ctl.load(ctlfname);

    StraightWG2D wg;
    clog << "Loading solution... ";
    wg.init( ctl ); // Initialize the simulation
    clog << "done\n";

    if ( useFFT )
    {
      arma::mat matrix;
      clog << "Evaluating field inside waveguide... ";
      wg.getFieldInsideWG( matrix );
      clog << "done\n";

      clog << "Computing FFT... ";
      arma::cx_mat fftField = arma::fft(matrix);
      matrix = abs(fftField);
      clog << "done\n";

      stringstream fname;
      fname << "data/fieldFFT" << ctl.getUID() << ".h5";
      matrix.save(fname.str().c_str(), arma::hdf5_binary);
      clog << "FFT written to " << fname.str() << endl;
    }
    else
    {
      double x0 = 0.0;
      double stepX = wg.transverseDiscretization().step;
      double stepZ = wg.longitudinalDiscretization().step;
      double width = wg.getWidth();
      double x = x0;

      vector<HarmInvRes> allRes;
      arma::mat fieldInside;

      // TODO: Now a copy of the real part of the solution is stored int the fieldInside matrix.
      //       This is not memory efficient. Optimize?
      wg.getFieldInsideWG( fieldInside );
      //while ( x < x0+width )
      for ( unsigned int k=0;k<fieldInside.n_rows;k++ )
      {
        clog << "Current x: " << x << " . Ends at: " << x0+width << endl;
        HarmInvRes results;
        vector<cdouble> res;
        vector<double> realRes;

        results.wcrd = x;

        // The harminv library requires a complex array. If real the complex part should be set to zero
        // Thus, copy the real field onto a complex array
        for ( unsigned int i=0;i<fieldInside.n_cols;i++ )
        {
          res.push_back(fieldInside(k,i));
          results.data.push_back(res[i].real());
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
          results.freq.push_back(harminv_get_freq( data, i)/stepZ);
          results.decay.push_back(harminv_get_decay( data, i)/stepZ);
          cdouble amplitudeCompx;
          amplitudeCompx = harminv_get_amplitude(data, i);
          results.amplitude.push_back(abs(amplitudeCompx));
          results.phase.push_back( arg(amplitudeCompx) );
        }
        x += stepX;
        allRes.push_back(results);
      }

      stringstream fname;
      fname << "data/harminvdata" << ctl.getUID() << ".h5";
      hid_t file_id = H5Fcreate(fname.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
      for ( unsigned int i=0;i<allRes.size();i++ )
      {
        stringstream freqname, decayname, ampname, phasename, dataname;
        freqname << "freq" << i;
        decayname << "decay" << i;
        ampname << "amplitude" << i;
        phasename << "phase" << i;
        dataname << "transformedData" << i;
        hsize_t dim = allRes[i].freq.size();
        unsigned int rank = 1;

        H5LTmake_dataset( file_id, freqname.str().c_str(), rank, &dim, H5T_NATIVE_DOUBLE, &allRes[i].freq[0]);
        H5LTmake_dataset( file_id, decayname.str().c_str(), rank, &dim, H5T_NATIVE_DOUBLE, &allRes[i].decay[0]);
        H5LTmake_dataset( file_id, ampname.str().c_str(), rank, &dim, H5T_NATIVE_DOUBLE, &allRes[i].amplitude[0]);
        H5LTmake_dataset( file_id, phasename.str().c_str(), rank, &dim, H5T_NATIVE_DOUBLE, &allRes[i].phase[0]);
        dim = allRes[i].data.size();
        H5LTmake_dataset( file_id, dataname.str().c_str(), rank, &dim, H5T_NATIVE_DOUBLE, &allRes[i].data[0]);

        // Set the attribute for the position
        H5LTset_attribute_double( file_id, freqname.str().c_str(), "x", &allRes[i].wcrd, 1);
        H5LTset_attribute_double( file_id, decayname.str().c_str(), "x", &allRes[i].wcrd, 1);
        H5LTset_attribute_double( file_id, ampname.str().c_str(), "x", &allRes[i].wcrd, 1);
        H5LTset_attribute_double( file_id, phasename.str().c_str(), "x", &allRes[i].wcrd, 1);
        H5LTset_attribute_double( file_id, dataname.str().c_str(), "x", &allRes[i].wcrd, 1);
      }
      H5Fclose(file_id);
      clog << "Data written to " << fname.str() << endl;
    }
  }
  catch ( exception &exc )
  {
    cout << exc.what() << endl;
    return 1;
  }
  return 0;
}
