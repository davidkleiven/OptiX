#include "refractiveIndex.hpp"
#include <stdexcept>
#include <fstream>
#include <cassert>
using namespace std;
const map<elementName, filename> RefractiveIndex::knownElements = {{"Ta", "MatProp/indexRefrTa.txt"},
                                    {"C2H6O", "MatProp/indexRefrC2H6O.txt"},
                                    {"SiO2", "MatProp/indexRefrSiO2.txt"},
                                    {"C2H6O2", "MatProp/indexRefrC2H6O2.txt"},
                                    {"Pb", "MatProp/indexRefrPb.txt"},
                                    {"Au", "MatProp/indexRefrAu.txt"},
                                    {"Vacuum", ""}};

void RefractiveIndex::load( const char* element )
{
  if ( knownElements.find(element) == knownElements.end() )
  {
    string msg("Unknown element ");
    msg += element;
    msg += "!";
    throw ( runtime_error(msg) );
  }
  string lookUp(element);
  if ( lookUp == "Vacuum")
  {
    isVacuum = true;
    return;
  }

  string fname(knownElements.at(lookUp));
  ifstream infile;
  infile.open( fname.c_str() );
  if ( !infile.good() )
  {
    throw (runtime_error("Could not open refractive index file!"));
  }

  // Read the first two lines
  string line;
  getline( infile, line);
  getline( infile, line);

  double newEnergy, newDelta, newBeta;
  while( infile >> newEnergy >> newDelta >> newBeta )
  {
    energy.push_back( newEnergy );
    delta.push_back( newDelta );
    beta.push_back( newBeta );
  }
  infile.close();
}

unsigned int RefractiveIndex::closestAbove( double energyInEv ) const
{
  for ( unsigned int i=0;i<energy.size();i++ )
  {
    if ( energyInEv < energy[i] )
    {
      return i;
    }
  }
  return energy.size();
}

double RefractiveIndex::getDelta( double energyInEv ) const
{
  if ( isVacuum ) return 0.0;

  unsigned int closest = closestAbove( energyInEv );
  if ( closest == energy.size() )
  {
    return delta.back();
  }
  else if ( closest == 0 )
  {
    return delta[0];
  }
  double weight = energyWeight( energyInEv, closest );
  return weight*delta[closest] + (1.0-weight)*delta[closest-1];
}

double RefractiveIndex::getBeta( double energyInEv ) const
{
  if ( isVacuum ) return 0.0;

  unsigned int closest = closestAbove( energyInEv );
  if ( closest == energy.size() )
  {
    return beta.back();
  }
  else if ( closest == 0 )
  {
    return beta[0];
  }
  double weight = energyWeight( energyInEv, closest);
  return weight*beta[closest] + (1.0-weight)*beta[closest-1];
}

double RefractiveIndex::energyWeight( double energyInEv, unsigned int closestIndxAbove  ) const
{
  double weight = (energyInEv-energy[closestIndxAbove-1])/( energy[closestIndxAbove]-energy[closestIndxAbove-1]);
  // DEBUG
  assert ( weight > 0.0 );
  assert ( weight < 1.0 );
  return weight;
}
