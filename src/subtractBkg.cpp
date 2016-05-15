#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

typedef vector<double> dvec;
int main(int argc, char** argv)
{
  if ( argc != 4 )
  {
    cout << "Usage: ./subtractBackground.cpp <ifile> <bkgfile> <ofile>\n";
    return 1;
  }

  string ifile(argv[1]);
  string bkgfile(argv[2]);
  string ofile(argv[3]);  

  dvec timeIn;
  dvec timeBkg;
  dvec realIn;
  dvec realBkg;
  dvec imagIn;
  dvec imagBkg;

  // Read data from infile
  string line;
  ifstream is(ifile.c_str());
  if ( !is.good() )
  {
    cout << "Problem when opening file " << ifile << endl;
    return 1;
  }

  while ( getline(is, line) && (line[0] == '#') ){};

  double time, realPart, imagPart;
  char comma;
  while ( is >> time  >> comma >> realPart >> comma >> imagPart )
  {
    timeIn.push_back(time);
    realIn.push_back(realPart);
    imagIn.push_back(imagPart);
  }
  is.close();

  // Read data from bkgfile
  is.open(bkgfile.c_str());
  
  if ( !is.good() )
  {
    cout << "Problem when opening file " << bkgfile << endl;
    return 1;
  }

  while ( getline(is, line) && (line[0] == '#') ){};

  while ( is >> time >> comma >> realPart >> comma >> imagPart )
  {
    timeBkg.push_back(time);
    realBkg.push_back(realPart);
    imagBkg.push_back(imagPart);
  }
  is.close();

  // Subtract off the background
  if ( timeIn.size() != timeBkg.size() )
  {
    cout << "The background time and the bkg time must have exactly the same timepoints!\n";
    return 1;
  }

  for ( unsigned int i=0;i<timeIn.size();i++)
  {
    realIn[i] -= realBkg[i];
    imagIn[i] -= imagBkg[i];
  }
  
  // Write the results to file
  ofstream os(ofile.c_str());
  if ( !os.good() )
  {
    cout << "Problems when opening file " << ofile << endl;
    return 1;
  }

  os << "# Field after subtracting of background\n";
  os << "# Infile: " << ifile << "\n";
  os << "# Background file " << bkgfile << "\n";
  os << "# Time, Ez.real, Ez.imag\n";
  for ( unsigned int i=0;i<timeIn.size();i++)
  {
    os << timeIn[i] << "," << realIn[i] << "," << imagIn[i] << "\n";
  }
  os.close();
  return 0;
}
