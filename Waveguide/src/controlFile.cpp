#include "controlFile.hpp"
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

using namespace std;

ControlFile::ControlFile( const char* name ):fname(name), base(new Json::Value), uid(0)
{
  srand(time(NULL));
  unsigned int uidMax = pow( 10, uidDigits );
  uid = rand()%uidMax;
  if ( uid < 100000 )
  {
    uid *= 10;
  }

  clog << "UID: " << uid << endl;

  (*base)["UID"] = uid;
  stringstream ss;
  ss << fname << uid;
  fname = ss.str();
}

ControlFile::~ControlFile()
{
  delete base;
}

void ControlFile::save() const
{
  string ofname = fname+".json";
  ofstream out(ofname);
  if ( !out.good() )
  {
    cout << "Could not open file " << ofname << endl;
    return;
  }

  Json::StyledWriter sw;
  out << sw.write(*base) << endl;
  out.close();
}
