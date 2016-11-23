#ifndef REFRACTIVE_INDEX_H
#define REFRACTIVE_INDEX_H
#include <string>
#include <map>
#include <vector>

typedef std::string filename;
typedef std::string elementName;
class RefractiveIndex
{
public:
  RefractiveIndex(){};
  void load( const char *element );
  double getDelta( double energyInEv ) const;
  double getBeta( double energyInEv ) const;
private:
  unsigned int closestAbove( double energyInEv ) const;
  double energyWeight( double energyInEv, unsigned int closestIndxAbove ) const;
  static const std::map<elementName, filename> knownElements;
  std::vector<double> energy;
  std::vector<double> delta;
  std::vector<double> beta;
};
#endif
