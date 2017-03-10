#include "refractiveIndexMaterial.hpp"
#include <json/reader.h>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <cassert>
#include <iostream>

using namespace std;

void SellmeierMaterial::load( const char* infname )
{
  fname = infname;

  Json::Value root;
  Json::Reader reader;

  ifstream infile;
  infile.open( fname.c_str() );
  if ( !infile.good() )
  {
    stringstream msg;
    msg << "Could not open file " << fname;
    throw( runtime_error( msg.str() ));
  }
  reader.parse( infile, root );
  infile.close();
  checkRequiredFields( root );

  for ( unsigned int i=0;i<root["lorentzians"].size();i++ )
  {
    Lorentzian lor;
    lor.preFactor = root["lorentzians"][i]["prefactor"].asDouble();
    lor.resonanceWavelength = root["lorentzians"][i]["resonance"].asDouble();
    lorentzians.push_back( lor );
  }
  epsInf = root["epsInf"].asDouble();
}

void SellmeierMaterial::checkRequiredFields( const Json::Value &root ) const
{
  if ( !root.isMember("epsInf") )
  {
    throw( runtime_error("No field named epsInf!") );
  }
  else if ( !root.isMember("lorentzians") )
  {
    throw( runtime_error("No field named lorentzians!") );
  }

  // Check all lorentzians
  for ( unsigned int i=0;i<root["lorentzians"].size();i++ )
  {
    if ( !root["lorentzians"][i].isMember("prefactor") )
    {
      stringstream msg;
      msg << "No field named prefactor in element " << i << " of the lorentzians!";
      throw( runtime_error(msg.str()));
    }
    else if ( !root["lorentzians"][i].isMember("resonance") )
    {
      stringstream msg;
      msg << "No field named resonance in element " << i << " of the lorentzians!";
      throw( runtime_error(msg.str()));
    }
  }
}

void SellmeierMaterial::getMEEPLorentzian( double lengthscaleInMicroMeter, unsigned int indx, double &sigma, double &resonnanceAngularFreq ) const
{
  assert( indx < lorentzians.size() );
  double PI = acos(-1.0);

  sigma = lorentzians[indx].preFactor;
  resonnanceAngularFreq = pow( 2.0*PI*lengthscaleInMicroMeter, 2 )/lorentzians[indx].resonanceWavelength;
}

double SellmeierMaterial::getEpsilon( double lambda ) const
{
  double eps = epsInf;
  for ( unsigned int i=0;i<lorentzians.size();i++ )
  {
    eps += lorentzians[i].preFactor*lambda*lambda/( lambda*lambda - lorentzians[i].resonanceWavelength );
  }
  return eps;
}

double SellmeierMaterial::getEpsilon( double lengthscale, double MEEPangFreq ) const
{
  const double PI = acos(-1.0);
  double lambda = 2.0*PI*lengthscale/MEEPangFreq;
  return getEpsilon( lambda );
}
