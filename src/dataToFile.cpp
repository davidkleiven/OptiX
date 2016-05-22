#include "dataToFile.h"
#include <stdexcept>
#include <sstream>
#include <fstream>

void monitorToFile(const std::string &fname, const dvec &time, const dvec &fieldVal, const meep::vec &pos)
{ 
  std::ofstream os(fname.c_str());
  if ( !os.good() )
  {
    std::stringstream ss;
    ss << "Error when opening file " << fname;
    throw (std::runtime_error(ss.str()));
  }
  os << "# Field monitored at positions\n";
  os << "# Position x=" << pos.x() << ",y=" << pos.y() << ",z="<<pos.z()<< std::endl;
  os << "# Time, Ez.real\n";
  for ( unsigned int i=0;i<time.size();i++ )
  {
    os << time[i] << "," << fieldVal[i] << "\n";
  }
  os.close(); 
}
