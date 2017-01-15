#include <pei/dialogBox.hpp>
#include "multipleCurvedWG.hpp"
#include "controlFile.hpp"
#include <map>
#include <string>

using namespace std;

int main( int argc, char** argv )
{
  map<string,double> params;
  params["wavelength"] = 0.1569;
  params["delta"] = 4.14E-5;
  params["beta"] = 4.45E-6;
  params["stepX"] = 0.5;
  params["stepZ"] = 80.0;
  params["downSamplingZ"] = 10.0;
  params["downSamplingX"] = 10;
  params["width"] = 100.0;
  params["nWaveguides"] = 3;
  params["Rmin"] = 0.01;
  params["maxIter"] = 1000;
  params["totalDeflectionAngle"] = 20.0;

  MultipleCurvedWG simulation;

  return 0;
}
