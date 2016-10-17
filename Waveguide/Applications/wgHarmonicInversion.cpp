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

using namespace std;
typedef complex<double> cdouble;

struct HarmInvRes
{
  vector<double> freq;
  vector<cdouble> amplitude;
  vector<double> decay;
  double wcrd;
};

int main( int argc, char** argv )
{
  string ctlfname("data/singleCurvedWG702004.json");
  unsigned int nFreq =300;
  double freqMin = 0.0;
  double freqMax = 0.1;

  ControlFile ctl;
  ctl.load(ctlfname);

  StraightWG2D wg;
  wg.init( ctl ); // Initialize the simulation

  double x0 = 0.0;
  double stepX = wg.transverseDiscretization().step;
  double width = wg.getWidth();

  double x = x0;

  vector<cdouble> res;
  while ( x < x0+width )
  {
    clog << "Current x: " << x << " . Ends at: " << x0+width << endl;
    HarmInvRes results;
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
    for ( unsigned int i=0;i<foundFreq;i++ )
    {
      results.freq.push_back(harminv_get_freq( data, i));
      results.decay.push_back(harminv_get_decay( data, i));
      cdouble amplitudeCompx;
      amplitudeCompx = harminv_get_amplitude(data, i);
      results.amplitude.push_back(amplitudeCompx);
    }
    x += stepX;
  }
  return 0;
}
