#include "controlFile.hpp"
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <jsoncpp/json/reader.h>

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

ControlFile::ControlFile(): base(new Json::Value), uid(0), fname(""){};

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

void ControlFile::load( const string &infname )
{
  ifstream infile(infname.c_str());
  if ( !infile.good() )
  {
    throw (runtime_error("Could not open control file"));
  }

  Json::Reader reader;
  reader.parse(infile, *base);
  uid = (*base)["UID"].asInt();

  // Remove ending which should be .json
  fname = infname.substr(0, infname.size()-4);
}
