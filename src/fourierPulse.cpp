#include <fstream>
#include <iostream>
#include <gsl/gsl_fft_real.h>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstring>
#define MIN_SAVE_VAL 0.0

using namespace std;
int main(int argc, char** argv)
{
  if ( argc != 2 )
  {
    cout << "Usage ./fourierPulse.out <pulsefile>\n";
    return 1;
  }

  // Read file
  vector<double> time;
  vector<double> field;
  string fname(argv[1]);
  ifstream infile(fname.c_str()); 
  if ( !infile.good() )
  {
    cout << "Problem when opening file " << fname << endl;
    return 1;
  }

  string line;
  while ( getline(infile, line) && (line[0] == '#') ){}
  stringstream firstline;
  firstline << line;
  double newT, newF;
  firstline >> newT >> newF;
  
  time.push_back(newT);
  field.push_back(newF);

  char comma;
  while ( infile >> newT >> comma >> newF )
  {
    time.push_back(newT);
    field.push_back(newF);
  }
  infile.close();

  // Compute fourier transform of the signal
  double dt = time[1] - time[0];
  
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
    field[field.size()/2] *= field[field.size()/2];
  }
  
  // Write results to file
  double df = 1.0/dt;
  unsigned int pos = fname.find(".");
  string ofname = fname.substr(0,pos);
  ofname += "Fourier.csv";
  ofstream os(ofname.c_str());
  os << "# Fourier transformed signal\n";
  os << "# Freq, Power\n";
  double freq = 0.0;
  for ( unsigned int i=0;i<field.size()/2;i++)
  {
    if ( field[i] > MIN_SAVE_VAL )
    {
      os << freq << "," << field[i] << "\n";
    }
    freq += df;
  }
  os.close();
  return 0;
}
  
      

  
  
  
