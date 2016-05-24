#include "readCSVdata.h"
#include <fstream>
#include <stdexcept>
#include <sstream>

using namespace std;
void ReadCSVData::read(const string &fname, unsigned int nColumn)
{
  nCol = nColumn;
  data.clear();
  ifstream infile(fname.c_str());

  if ( !infile.good() )
  {
    stringstream msg;
    msg << "Could not open file " << fname;
    throw ( runtime_error(msg.str()) );
  }
  
  // Read header
  string line;
  while ( getline(infile, line) )
  {
    if ( line[0] == '#' )
    {
      continue;
    }
    stringstream lineRead;
    lineRead << line;
    double newVal;
    lineRead >> newVal;
    data.push_back(newVal);
    for ( unsigned int i=0;i<nCol-1;i++)
    {
      char comma;
      lineRead >> comma >> newVal;
      data.push_back(newVal);
    }
  }
  infile.close();
}
      
double ReadCSVData::get(unsigned int row, unsigned int col) const
{
  unsigned int indx = row*nCol + col;
  return data[indx];
}
  
