#include <fstream>
#include <iostream>
#include <gsl/gsl_fft_real.h>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstring>
#include "readCSVdata.h"
#define MIN_SAVE_VAL 1E-4

using namespace std;
int main(int argc, char** argv)
{
  if ( argc != 2 )
  {
    cout << "Usage ./fourierPulse.out <pulsefile>\n";
    return 1;
  }

  // Read file
  string fname(argv[1]);
  ReadCSVData reader;
  reader.read(fname, 2);
  // Compute sum
  double fieldSum;
  for ( unsigned int i=0;i<reader.numPoints();i++)
  {
    fieldSum += reader.get(i,1)*reader.get(i,1);
  }

  // Compute fourier transform of the signal
  double dt = reader.get(1,0) - reader.get(0,0);

  vector<double> field;
  // Store copy of field in array for FFT
  for ( unsigned int i=0;i<reader.numPoints();i++ )
  {
    field.push_back(reader.get(i,1));
  }
  
  gsl_fft_real_wavetable *realTab;
  gsl_fft_real_workspace *work;
  
  work = gsl_fft_real_workspace_alloc(field.size());
  realTab = gsl_fft_real_wavetable_alloc(field.size());
  gsl_fft_real_transform(&field[0], 1, field.size(), realTab, work);

  gsl_fft_real_wavetable_free(realTab);
  gsl_fft_real_workspace_free(work);

  if ( field.size()%2 == 1 )
  {
    // Odd length
    field[0] *= field[0];
    for ( unsigned int i=1;i<field.size()/2;i++ )
    {
      field[i] = pow( field[2*i-1], 2 ) + pow( field[2*i], 2 );
    }
  }
  else
  {
    // Even length
    field[0] *= field[0];
    for ( unsigned int i=1;i<field.size()/2-1;i++)
    {
      field[i] = pow( field[2*i-1], 2 ) + pow( field[2*i], 2 );
    }
    field[field.size()/2-1] *= field[field.size()/2];
  }

  // Compute sum of spectrum
  double fieldFFTSum = 0.0;
  for ( unsigned int i=0;i<field.size()/2;i++)
  {
    fieldFFTSum += field[i];
  }
  fieldFFTSum = 2.0*fieldFFTSum - field[0];
  
  // Check Parseval
  if ( abs(fieldFFTSum/static_cast<double>(field.size()) - fieldSum) > 1E-3 )
  {
    cout << "Warning! Parseval identity is not fulfilled\n";
    cout << "Sum of field: " << fieldSum << endl;
    cout << "Sum of power spectrum (/N): " << fieldFFTSum/static_cast<double>(field.size()) << endl;
  }
  
  // Write results to file
  double df = 1.0/static_cast<double>(field.size());
  unsigned int pos = fname.find(".");
  string ofname = fname.substr(0,pos);
  ofname += "Fourier.csv";
  ofstream os(ofname.c_str());
  os << "# Fourier transformed signal\n";
  os << "# Freq, Power\n";
  double freq = 0.0;
  for ( unsigned int i=0;i<field.size()/2;i++)
  {
    double fieldSave = field[i];//static_cast<double>(field.size());
    if ( fieldSave > MIN_SAVE_VAL )
    {
      os << freq << "," << fieldSave << "\n";
    }
    freq += df;
  }
  os.close();
  return 0;
} 
