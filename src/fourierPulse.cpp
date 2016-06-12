#include <fstream>
#include <iostream>
#include <gsl/gsl_fft_complex.h>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstring>
#include "readCSVdata.h"
#include <complex>
#include <stdexcept>
#include <jsoncpp/json/reader.h>
#include <jsoncpp/json/writer.h>
#define MIN_RELATIVE_SAVE_VAL 0.05
//#define OUTPUT_SUBTRACTED
/**
* This file takes three command line arguments
* 1) name of file with scatterer
* 2) name of file without scatterer
* 3) Incident angle in degrees of peak frequency
*
* Example:
* ./fourierPulse.out fileWithSc.csv fileWithoutSc.csv 20
*
* The format of the csv files are assumed to be
* <time> <field at transmission point> <field at reflection point>
*/
const double PI = acos(-1.0);
using namespace std;
int main(int argc, char** argv)
{
  string id("[FFT field] ");
  if ( argc != 4 )
  {
    cout << id << "Usage ./fourierPulse.out <infile> <bkgfile> <incident angle of peak>\n";
    return 1;
  }

  // Read file
  string fname(argv[1]);
  string bkgfname(argv[2]);
  stringstream ss;  
  ss << argv[3];
  double angle;
  ss >> angle;

  Json::Reader runReader;
  Json::Reader bkgReader;

  ifstream bkgstream(bkgfname.c_str());
  if ( !bkgstream.good() )
  {
    cerr << id << "Could not open file " << bkgfname << endl;
    return 1;
  }
  Json::Value bkg;
  bkgReader.parse( bkgstream, bkg );
  bkgstream.close();

  ifstream runStream( fname.c_str() );
  if ( !runStream.good() )
  {
    cerr << id << "Could not open file " << fname << endl;
    return 1;
  }
  Json::Value run;
  runReader.parse( runStream, run );
  runStream.close();
  
  vector<double> fieldRefl;
  unsigned int nPointsBkg = bkg["timepoints"].size();
  unsigned int nPointsRun = run["timepoints"].size();
  if ( nPointsBkg != nPointsRun )
  {
    cout << "Different number of points in infile and bkg file\n",  
    cout << "Number of points in infile: " << nPointsBkg << endl;
    cout << "Number of points in bkgfile: " << nPointsRun << endl;
    return 1;
  }

  // Subtract off difference
  for ( unsigned int i=0;i<nPointsBkg;i++ )
  {
    run["reflected"]["real"][i] = run["reflected"]["real"][i].asDouble() - bkg["reflected"]["real"][i].asDouble();
    run["reflected"]["imag"][i] = run["reflected"]["imag"][i].asDouble() - bkg["reflected"]["imag"][i].asDouble();
  }

  // Compute fourier transform of the signal
  double dt = run["timepoints"][1].asDouble() - run["timepoints"][0].asDouble();

  // Store arrays as required by the FFT routines
  vector<double> reflectedBkg;
  vector<double> transmittedBkg;
  vector<double> reflectedRun;
  vector<double> transmittedRun;
  for ( unsigned int i=0;i<nPointsBkg;i++ )
  {
    reflectedBkg.push_back( bkg["reflected"]["real"][i].asDouble() );
    reflectedBkg.push_back( bkg["reflected"]["imag"][i].asDouble() );
    transmittedBkg.push_back( bkg["transmitted"]["real"][i].asDouble() );
    transmittedBkg.push_back( bkg["transmitted"]["imag"][i].asDouble() );

    reflectedRun.push_back( run["reflected"]["real"][i].asDouble() );
    reflectedRun.push_back( run["reflected"]["imag"][i].asDouble() );
    transmittedRun.push_back( run["transmitted"]["real"][i].asDouble() );
    transmittedRun.push_back( run["transmitted"]["imag"][i].asDouble() );
  }
  
  // Coompute Fourier transforms
  gsl_fft_complex_wavetable *cTab;
  gsl_fft_complex_workspace *work;
  
  work = gsl_fft_complex_workspace_alloc(nPointsBkg);
  cTab = gsl_fft_complex_wavetable_alloc(nPointsBkg);
  gsl_fft_complex_forward(&reflectedBkg[0], 1, nPointsBkg, cTab, work);
  gsl_fft_complex_forward(&transmittedBkg[0], 1, nPointsBkg, cTab, work); 
  gsl_fft_complex_forward(&reflectedRun[0], 1, nPointsBkg, cTab, work);
  gsl_fft_complex_forward(&transmittedRun[0], 1, nPointsBkg, cTab, work);

  gsl_fft_complex_wavetable_free(cTab);
  gsl_fft_complex_workspace_free(work);

  // Find position of maximum use the one from the transmitted signal
  unsigned int currentMaxPos = 1;
  double currentMax = -1.0;
  for ( unsigned int i=0;i<nPointsBkg;i++)
  {
    double normTrans = sqrt( transmittedRun[2*i]*transmittedRun[2*i] + transmittedRun[2*i+1]*transmittedRun[2*i+1] );
    if ( normTrans > currentMax )
    {
      currentMaxPos = i;
      currentMax = normTrans;
    }
  }

  Json::Value transCoeffReal(Json::arrayValue);
  Json::Value transCoeffImag(Json::arrayValue);
  Json::Value reflCoeffReal(Json::arrayValue);
  Json::Value reflCoeffImag(Json::arrayValue);
  Json::Value angleArray(Json::arrayValue);
  double df = 1.0/(dt*static_cast<double>(nPointsBkg));
  double frequencyAtMax = static_cast<double>(currentMaxPos)*df;
  double estimateOfMaxTransValue = sqrt(2.0)*currentMax;
  for ( unsigned int i=0;i<nPointsBkg;i++ )
  {
    complex<double> refl(reflectedRun[2*i], reflectedRun[2*i+1]);
    complex<double> bkgVal(reflectedBkg[2*i], reflectedBkg[2*i+1]);
    complex<double> refCoeff = refl/bkgVal;
    
    complex<double> trans(transmittedRun[2*i], transmittedRun[2*i+1]);
    bkgVal = (transmittedBkg[2*i], transmittedBkg[2*i+1]);
    complex<double> transCoeff = trans/bkgVal;

    double freq = static_cast<double>(i)/(dt*static_cast<double>(nPointsBkg));
    double angleArgument = frequencyAtMax*sin( angle*PI/180.0 )/freq;
  
    if ( abs(angleArgument) > 1.0 )
    {
      continue;
    }
    double currentAngle = asin( angleArgument )*180.0/PI;
    // Compute norm of reflected and transmitted fields and save only the significant
    double reflNorm = abs( refCoeff );
    double transNorm = abs( transCoeff );
    if ( (reflNorm > MIN_RELATIVE_SAVE_VAL*estimateOfMaxTransValue) || \
       ( transNorm > MIN_RELATIVE_SAVE_VAL*estimateOfMaxTransValue) )
    {
      reflCoeffReal.append( real(refCoeff) );
      reflCoeffImag.append( imag(refCoeff) );
    
      transCoeffReal.append( real( transCoeff ) );
      transCoeffImag.append( imag( transCoeff ) );
      angleArray.append( currentAngle );
    }
  }
 
  // Write results to file
  Json::Value base;
  base["geometry"] = run["geometry"];
  base["reflection"]["real"] = reflCoeffReal;
  base["reflection"]["imag"] = reflCoeffImag;
  base["transmition"]["real"] = transCoeffReal;
  base["transmition"]["imag"] = transCoeffImag;
  base["angle"] = angleArray;
  Json::FastWriter fw;

  unsigned int pos = fname.find(".");
  string ofname = fname.substr(0,pos);
  ofname += "Fourier.csv";
  ofstream os(ofname.c_str());

  if ( !os.good() )
  {
    cerr << id << "Could not open file " << ofname << endl;
    return 1;
  }
  os << fw.write( base ) << endl;
  os.close();
  return 0;
} 
