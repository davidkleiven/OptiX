#include <fstream>
#include <iostream>
#include <gsl/gsl_fft_real.h>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstring>
#include "readCSVdata.h"
#include <complex>
#include <stdexcept>
#define MIN_SAVE_VAL 1E-4
//#define OUTPUT_SUBTRACTED

const double PI = acos(-1.0);
using namespace std;
int main(int argc, char** argv)
{
  if ( argc != 4 )
  {
    cout << "Usage ./fourierPulse.out <infile> <bkgfile> <incident angle of peak>\n";
    return 1;
  }

  // Read file
  string fname(argv[1]);
  string bkgfname(argv[2]);
  stringstream ss;  
  ss << argv[3];
  double angle;
  ss >> angle;

  ReadCSVData reader;
  ReadCSVData bkgreader;
  try
  {
    reader.read(fname, 3);
    bkgreader.read(bkgfname, 3);
  }
  catch(runtime_error &exc)
  {
    cerr << exc.what() << endl;
    return 1;
  }
  catch(...)
  {
    cerr << "Unknown exception...\n";
    return 1;
  }

  vector<double> fieldRefl;
  if ( bkgreader.numPoints() != reader.numPoints() )
  {
    cout << "Different number of points in infile and bkg file\n",  
    cout << "Number of points in infile: " << reader.numPoints() << endl;
    cout << "Number of points in bkgfile: " << bkgreader.numPoints() << endl;
    return 1;
  }

  // Subtract off difference
  for ( unsigned int i=0;i<bkgreader.numPoints();i++ )
  {
    fieldRefl.push_back( reader.get(0,2) - bkgreader.get(0,2) );
  }

  // Compute sum
  double fieldSumTrans = 0.0;
  double fieldSumRefl = 0.0;
  for ( unsigned int i=0;i<reader.numPoints();i++)
  {
    fieldSumTrans += reader.get(i,1)*reader.get(i,1);
    fieldSumRefl += fieldRefl[i]*fieldRefl[i];
  }

  // Compute fourier transform of the signal
  double dt = reader.get(1,0) - reader.get(0,0);

  vector<double> fieldTrans;
  vector<double> bkgRefl;
  vector<double> bkgTrans;
  // Store copy of field in array for FFT
  for ( unsigned int i=0;i<reader.numPoints();i++ )
  {
    fieldTrans.push_back(reader.get(i,1));
    bkgTrans.push_back(bkgreader.get(i,1));
    bkgRefl.push_back(bkgreader.get(i,2));
  }
  
  // Coompute Fourier transforms
  gsl_fft_real_wavetable *realTab;
  gsl_fft_real_workspace *work;
  
  work = gsl_fft_real_workspace_alloc(fieldTrans.size());
  realTab = gsl_fft_real_wavetable_alloc(fieldTrans.size());
  gsl_fft_real_transform(&fieldTrans[0], 1, fieldTrans.size(), realTab, work);
  gsl_fft_real_transform(&fieldRefl[0], 1, fieldRefl.size(), realTab, work); 
  gsl_fft_real_transform(&bkgTrans[0], 1, bkgTrans.size(), realTab, work);
  gsl_fft_real_transform(&bkgRefl[0], 1, bkgRefl.size(), realTab, work);

  gsl_fft_real_wavetable_free(realTab);
  gsl_fft_real_workspace_free(work);

  vector< complex<double> > reflection;
  vector< complex<double> > transmission;
  if ( fieldTrans.size()%2 == 1 )
  {
    // Odd length
    fieldTrans[0] *= fieldTrans[0];
    fieldRefl[0] *= fieldRefl[0];
    reflection.push_back( (fieldRefl[0]/bkgRefl[0],0.0) );
    for ( unsigned int i=1;i<fieldTrans.size()/2;i++ )
    {
      fieldTrans[i] = pow( fieldTrans[2*i-1], 2 ) + pow( fieldTrans[2*i], 2 );
    }
    for ( unsigned int i=1;i<fieldRefl.size()/2;i++ )
    {
      fieldRefl[i] = pow( fieldRefl[2*i-1], 2 ) + pow( fieldRefl[2*i], 2 );
    }
  }
  else
  {
    // Even length
    fieldTrans[0] *= fieldTrans[0];
    fieldRefl[0] *= fieldRefl[0];
    for ( unsigned int i=1;i<fieldTrans.size()/2-1;i++)
    {
      fieldTrans[i] = pow( fieldTrans[2*i-1], 2 ) + pow( fieldTrans[2*i], 2 );
    }
    fieldTrans[fieldTrans.size()/2-1] *= fieldTrans[fieldTrans.size()/2];
    for ( unsigned int i=1;i<fieldRefl.size()/2-1;i++)
    {
      fieldRefl[i] = pow( fieldRefl[2*i-1], 2 ) + pow( fieldRefl[2*i], 2 );
    }
    fieldRefl[fieldRefl.size()/2-1] *= fieldRefl[fieldRefl.size()/2];
  }

  // Compute sum of spectrum
  double fieldTransFFTSum = 0.0;
  double fieldReflFFTSum = 0.0;
  for ( unsigned int i=0;i<fieldTrans.size()/2;i++)
  {
    fieldTransFFTSum += fieldTrans[i];
    fieldReflFFTSum += fieldRefl[i];
  }
  fieldTransFFTSum = 2.0*fieldTransFFTSum - fieldTrans[0];
  fieldReflFFTSum = 2.0*fieldReflFFTSum - fieldRefl[0];
  
  // Check Parseval
  if ( abs(fieldTransFFTSum/static_cast<double>(fieldTrans.size()) - fieldSumTrans) > 1E-3 )
  {
    cout << "Warning! Parseval identity for transmitted signal is not fulfilled\n";
    cout << "Sum of field: " << fieldSumTrans << endl;
    cout << "Sum of power spectrum (/N): " << fieldTransFFTSum/static_cast<double>(fieldTrans.size()) << endl;
  }
  if ( abs(fieldReflFFTSum/static_cast<double>(fieldRefl.size()) - fieldSumRefl) > 1E-3 )
  {
    cout << "Warning! Parseval identity for reflected signal is not fulfilled\n";
    cout << "Sum of field: " << fieldSumRefl << endl;
    cout << "Sum of power spectrum (/N): " << fieldReflFFTSum/static_cast<double>(fieldRefl.size()) << endl;
  }

  // Find position of maximum use the one from the transmitted signal
  unsigned int currentMaxPos = 1;
  double currentMax = -1.0;
  for ( unsigned int i=1;i<fieldTrans.size()/2;i++)
  {
    if ( fieldTrans[i] > currentMax )
    {
      currentMaxPos = i;
      currentMax = fieldTrans[i];
    }
  }
 
  // Write results to file
  double df = 1.0/(dt*static_cast<double>(fieldTrans.size()));
  double frequencyAtMax = static_cast<double>(currentMaxPos)*df;
  unsigned int pos = fname.find(".");
  string ofname = fname.substr(0,pos);
  ofname += "Fourier.csv";
  ofstream os(ofname.c_str());
  os << "# Fourier transformed signal\n";
  os << "# Freq, Transmitted, Reflected\n";
  double freq = 0.0;
  for ( unsigned int i=0;i<fieldTrans.size()/2;i++ )
  {
    double freq = static_cast<double>(i)/(dt*static_cast<double>(fieldTrans.size()));
    double angleArgument = 0.0;
    if ( freq > 0.0 )
    {
      angleArgument = frequencyAtMax*sin( angle*PI/180.0 )/freq;
    }
    else
    {
      continue;
    }
  
    if ( abs(angleArgument) > 1.0 )
    {
      continue;
    }
    double currentAngle = asin( angleArgument )*180.0/PI;
    double fieldSaveTrans = fieldTrans[i];//static_cast<double>(field.size());
    double fieldSaveRefl = fieldRefl[i];
    //if ( (fieldSaveTrans > MIN_SAVE_VAL) || (fieldSaveRefl > MIN_SAVE_VAL) )
    {
      os << freq << "," << currentAngle << "," << fieldSaveTrans << "," << fieldSaveRefl << "\n";
    }
  }
  os.close();
  return 0;
} 
