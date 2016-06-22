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
#include <gsl/gsl_fft_real.h>

#define MIN_RELATIVE_SAVE_VAL 1E-1
#define VERBOSE
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

    transmittedBkg.push_back( bkg["transmitted"]["real"][i].asDouble() );

    reflectedRun.push_back( run["reflected"]["real"][i].asDouble() );

    transmittedRun.push_back( run["transmitted"]["real"][i].asDouble() );
  }
  
  // Compute Fourier transforms
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

  // Find position of maximum use the one from the transmitted signal
  unsigned int currentMaxPos = 1;
  double currentMax = -1.0;
  unsigned int endIteration = nPointsBkg/2;
  for (unsigned int i=0;i<endIteration;i++)
  {
    double normTrans = normRealOddFFT(transmittedRun, i);
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
  Json::Value angleArrayT(Json::arrayValue);
  Json::Value freqArrayT(Json::arrayValue);
  Json::Value angleArrayR(Json::arrayValue);
  Json::Value freqArrayR(Json::arrayValue);
  double df = 1.0/(dt*static_cast<double>(nPointsBkg));
  double frequencyAtMax = static_cast<double>(currentMaxPos)*df;
  double estimateOfMaxTransValue = sqrt(2.0)*currentMax;
  double estimateOfMaxReflValue = sqrt(2.0)*normRealOddFFT(reflectedRun, currentMaxPos);

  #ifdef VERBOSE
    bool hasPrintedHFieldMsg = false;
    cout << "Field component: " << run["FieldComponent"].asString() << endl;
  #endif
  for ( unsigned int i=0;i<endIteration;i++ )
  {
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
    complex<double> refCoeff = refl/bkgVal;
    double refCoeffAngle = arg( refCoeff );

    // TODO: Handling of transmission coefficient is currently wrong
    if ( i > 0) 
    {
      realRun = transmittedRun[2*i-1];
      imagRun = transmittedRun[2*i];
      realBkg = transmittedBkg[2*i-1];
      imagBkg = transmittedBkg[2*i];
    }
    else
    {
      realRun = transmittedRun[0];
      imagRun = 0.0;  
      realBkg = transmittedBkg[0];
      imagBkg = 0.0;
    }
    complex<double> trans(realRun, imagRun);
    bkgVal.real(realBkg);
    bkgVal.imag(imagBkg);
      
    complex<double> transCoeff = trans/bkgVal;
    if ( run["FieldComponent"].asString() == "Hz" )
    {
      #ifdef VERBOSE
        if (!hasPrintedHFieldMsg)
        {
          cout << "Using H-field... Epsilon="<<run["geometry"]["EpsilonHigh"].asDouble() << endl;
          hasPrintedHFieldMsg = true;
        }
      #endif
      transCoeff /= sqrt( run["geometry"]["EpsilonHigh"].asDouble() );
    }
    double transCoeffAngle = arg( transCoeff );

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
    if (reflNorm > MIN_RELATIVE_SAVE_VAL*estimateOfMaxReflValue) 
    {
      reflCoeffNorm.append( abs(refCoeff) );
      reflCoeffPhase.append( refCoeffAngle );
    
      angleArrayR.append( currentAngle );
      freqArrayR.append( freq );
    }
    if ( transNorm > MIN_RELATIVE_SAVE_VAL*estimateOfMaxTransValue) 
    {
      angleArrayT.append( currentAngle );
      freqArrayT.append( freq );
      transCoeffNorm.append( abs(transCoeff)  );
      transCoeffPhase.append( transCoeffAngle );
    }
  }
 
  // Write results to file
  Json::Value base;
  base["geometry"] = run["geometry"];
  base["reflected"]["norm"] = reflCoeffNorm;
  base["reflected"]["phase"] = reflCoeffPhase;
  base["reflected"]["frequency"] = freqArrayR;
  base["reflected"]["angle"] = angleArrayR;
  base["reflected"]["position"] = run["reflected"]["position"];
  base["transmitted"]["norm"] = transCoeffNorm;
  base["transmitted"]["phase"] = transCoeffPhase;
  base["transmitted"]["frequency"] = freqArrayT;
  base["transmitted"]["angle"] = angleArrayT;
  base["transmitted"]["position"] = run["transmitted"]["position"];
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
