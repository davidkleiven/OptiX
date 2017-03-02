#include "paraxialSource.hpp"
#include <cmath>
#include <stdexcept>
#include <iostream>

using namespace std;

const double PI = acos(-1.0);

void ParaxialSource::setWavelength( double lambda )
{
  wavenumber = 2.0*PI/lambda;
}

double ParaxialSource::getWavelength() const
{
  return 2.0*PI/wavenumber;
}

double ParaxialSource::getWavenumber() const
{
  return wavenumber;
}

void ParaxialSource::info( Json::Value &obj ) const
{
  obj["name"] = name;
  obj["wavelength"] = getWavelength();
  obj["amplitude"] = amplitude;
  obj["z0"] = z0;
}

cdouble ParaxialSource::get( double x, double z ) const
{
  throw ( runtime_error("2D version of the ParaxialSorce get function is not implemented!") );
}

cdouble ParaxialSource::get( double x, double y, double z ) const
{
  throw( runtime_error("3D version of the ParaxialSource get function is not implemented!") );
}
