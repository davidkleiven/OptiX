#include <fstream>
#include <iostream>
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
//#define USE_COMPLEX_FIELD
#ifdef USE_COMPLEX_FIELD
  #include <gsl/gsl_fft_complex.h>
#else
  #include <gsl/gsl_fft_real.h>
#endif
#define MIN_RELATIVE_SAVE_VAL 1E-2
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

double normRealOddFFT(const vector<double> &fftOut, unsigned int indx)
{
  if ( indx == 0 )
  {
    return fftOut[0]*fftOut[0];
  }
  return sqrt( pow(fftOut[2*indx-1],2) + pow(fftOut[2*indx],2));
}

double normComplexFFT(const vector<double> &fftOut, unsigned int indx)
{
  return sqrt(pow(fftOut[2*indx],2) + pow(fftOut[2*indx+1],2));
}

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
  unsigned int nPointsBkg = bkg["time"].size();
  unsigned int nPointsRun = run["time"].size();
  if ( nPointsBkg != nPointsRun )
  {
    cout << "Different number of points in infile and bkg file\n",  
    cout << "Number of points in infile: " << nPointsBkg << endl;
    cout << "Number of points in bkgfile: " << nPointsRun << endl;
    return 1;
  }

  #ifndef USE_COMPLEX_FIELD
    // If real fields are used, make sure that the number of values are odd. 
    // Real FFT are different depending on if the number is even or odd.
    if ( nPointsBkg%2 == 0 )
    {
      nPointsBkg -= 1;
    }
  #endif

  // Subtract off difference
  for ( unsigned int i=0;i<nPointsBkg;i++ )
  {
    run["reflected"]["real"][i] = run["reflected"]["real"][i].asDouble() - bkg["reflected"]["real"][i].asDouble();
    #ifdef USE_COMPLEX_FIELD
      run["reflected"]["imag"][i] = run["reflected"]["imag"][i].asDouble() - bkg["reflected"]["imag"][i].asDouble();
    #endif
  }

  // Compute fourier transform of the signal
  double dt = run["time"][1].asDouble() - run["time"][0].asDouble();

  // Store arrays as required by the FFT routines
  vector<double> reflectedBkg;
  vector<double> transmittedBkg;
  vector<double> reflectedRun;
  vector<double> transmittedRun;
  for ( unsigned int i=0;i<nPointsBkg;i++ )
  {
    reflectedBkg.push_back( bkg["reflected"]["real"][i].asDouble() );
    #ifdef USE_COMPLEX_FIELD
      reflectedBkg.push_back( bkg["reflected"]["imag"][i].asDouble() );
    #endif
    transmittedBkg.push_back( bkg["transmitted"]["real"][i].asDouble() );
    #ifdef USE_COMPLEX_FIELD
      transmittedBkg.push_back( bkg["transmitted"]["imag"][i].asDouble() );
    #endif

    reflectedRun.push_back( run["reflected"]["real"][i].asDouble() );
    #ifdef USE_COMPLEX_FIELD
      reflectedRun.push_back( run["reflected"]["imag"][i].asDouble() );
    #endif
    transmittedRun.push_back( run["transmitted"]["real"][i].asDouble() );
    #ifdef USE_COMPLEX_FIELD
      transmittedRun.push_back( run["transmitted"]["imag"][i].asDouble() );
    #endif
  }
  
  // Coompute Fourier transforms
  #ifdef USE_COMPLEX_FIELD
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
  #else
    gsl_fft_real_wavetable *rTab;
    gsl_fft_real_workspace *work;
    
    work = gsl_fft_real_workspace_alloc(nPointsBkg);
    rTab = gsl_fft_real_wavetable_alloc(nPointsBkg);
    
    gsl_fft_real_transform(&reflectedBkg[0], 1, nPointsBkg, rTab, work);
    gsl_fft_real_transform(&transmittedBkg[0], 1, nPointsBkg, rTab, work); 
    gsl_fft_real_transform(&reflectedRun[0], 1, nPointsBkg, rTab, work);
    gsl_fft_real_transform(&transmittedRun[0], 1, nPointsBkg, rTab, work);

    gsl_fft_real_wavetable_free(rTab);
    gsl_fft_real_workspace_free(work);
  #endif

  // Find position of maximum use the one from the transmitted signal
  unsigned int currentMaxPos = 1;
  double currentMax = -1.0;
  #ifdef USE_COMPLEX_FIELD
    unsigned int endIteration = nPointsBkg;
  #else
    unsigned int endIteration = nPointsBkg/2;
  #endif
  for (unsigned int i=0;i<endIteration;i++)
  {
    #ifdef USE_COMPLEX_FIELD
      double normTrans = normComplexFFT(transmittedRun, i);
    #else
      double normTrans = normRealOddFFT(transmittedRun, i);
    #endif
    if ( normTrans > currentMax )
    {
      currentMaxPos = i;
      currentMax = normTrans;
    }
  }

  Json::Value transCoeffNorm(Json::arrayValue);
  Json::Value transCoeffPhase(Json::arrayValue);
  Json::Value reflCoeffNorm(Json::arrayValue);
  Json::Value reflCoeffPhase(Json::arrayValue);
  Json::Value angleArray(Json::arrayValue);
  Json::Value freqArray(Json::arrayValue);
  double df = 1.0/(dt*static_cast<double>(nPointsBkg));
  double frequencyAtMax = static_cast<double>(currentMaxPos)*df;
  double estimateOfMaxTransValue = sqrt(2.0)*currentMax;
  for ( unsigned int i=0;i<endIteration;i++ )
  {
    #ifdef USE_COMPLEX_FIELD
      complex<double> refl(reflectedRun[2*i], reflectedRun[2*i+1]);
      complex<double> bkgVal(reflectedBkg[2*i], reflectedBkg[2*i+1]);
    #else 
      double realRun, imagRun, realBkg, imagBkg;
      if ( i > 0) 
      {
        realRun = reflectedRun[2*i-1];
        imagRun = reflectedRun[2*i];
        realBkg = reflectedBkg[2*i-1];
        imagBkg = reflectedBkg[2*i];
      }
      else
      {
        realRun = reflectedRun[0];
        imagRun = 0.0;  
        realBkg = reflectedBkg[0];
        imagBkg = 0.0;
      }
      complex<double> refl(realRun, imagRun);
      complex<double> bkgVal(realBkg, imagBkg);
    #endif
    complex<double> refCoeff = refl/bkgVal;
    double refCoeffAngle = atan( imag(refCoeff)/real(refCoeff) );
    if (( imag(refCoeff) < 0.0 ) && ( real(refCoeff) < 0.0 ) )
    {
      // In third quadrant
      refCoeffAngle -= PI;
    }
    else if ( ( imag(refCoeff) > 0.0 ) && ( real(refCoeff) < 0.0 ))
    {
      // In second quadrant
      refCoeffAngle += PI;
    }
    
    // TODO: Handling of transmission coefficient is currently wrong
    complex<double> trans(transmittedRun[2*i], transmittedRun[2*i+1]);
    bkgVal = (transmittedBkg[2*i], transmittedBkg[2*i+1]);
    double transCoeff = abs( trans/bkgVal );
    double transCoeffAngle = atan( imag(trans)/real(trans) );

    double freq = static_cast<double>(i)/(dt*static_cast<double>(nPointsBkg));
    double angleArgument = frequencyAtMax*sin( angle*PI/180.0 )/freq;
  
    if ( abs(angleArgument) > 1.0 )
    {
      continue;
    }
    double currentAngle = asin( angleArgument )*180.0/PI;

    // Compute norm of reflected and transmitted fields and save only the significant
    double reflNorm = abs( refl );
    double transNorm = abs( trans );
    if ( (reflNorm > MIN_RELATIVE_SAVE_VAL*estimateOfMaxTransValue) || \
       ( transNorm > MIN_RELATIVE_SAVE_VAL*estimateOfMaxTransValue) )
    {
      reflCoeffNorm.append( abs(refCoeff) );
      reflCoeffPhase.append( refCoeffAngle );
    
      transCoeffNorm.append( transCoeff  );
      transCoeffPhase.append( transCoeffAngle );
      angleArray.append( currentAngle );
      freqArray.append( freq );
    }
  }
 
  // Write results to file
  Json::Value base;
  base["geometry"] = run["geometry"];
  base["reflection"]["norm"] = reflCoeffNorm;
  base["reflection"]["phase"] = reflCoeffPhase;
  base["transmition"]["norm"] = transCoeffNorm;
  base["transmition"]["phase"] = transCoeffPhase;
  base["angle"] = angleArray;
  base["frequency"] = freqArray;
  Json::FastWriter fw;

  unsigned int pos = fname.find(".");
  string ofname = fname.substr(0,pos);
  ofname += "Fourier.json";
  ofstream os(ofname.c_str());

  if ( !os.good() )
  {
    cerr << id << "Could not open file " << ofname << endl;
    return 1;
  }
  os << fw.write( base ) << endl;
  os.close();
  cout << id << "Coefficients written to " << ofname << endl;
  return 0;
} 
