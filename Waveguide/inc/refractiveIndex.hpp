#ifndef REFRACTIVE_INDEX_H
#define REFRACTIVE_INDEX_H
#include <string>
#include <map>
#include <vector>

typedef std::string filename;
typedef std::string elementName;

/** Class for handling energy dependent refractive indices */
class RefractiveIndex
{
public:
  RefractiveIndex(){};

  /** Load data from file for a given element */
  void load( const char *element );

  /** Get delta at the given energy in eV */
  double getDelta( double energyInEv ) const;

  /** Get beta at the given energy in eV */
  double getBeta( double energyInEv ) const;
private:
  /** Get closest index above the given energy in eV */
  unsigned int closestAbove( double energyInEv ) const;

  /** Get the weighting factor for linear interpolation */
  double energyWeight( double energyInEv, unsigned int closestIndxAbove ) const;
  static const std::map<elementName, filename> knownElements;
  std::vector<double> energy;
  std::vector<double> delta;
  std::vector<double> beta;
};
#endif
