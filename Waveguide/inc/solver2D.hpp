#ifndef SOLVER2D_H
#define SOLVER2D_H
#include<string>

class WaveGuideFDSimulation;

class Solver2D
{
public:
  Solver2D( const char* name ):name(name){};
  virtual ~Solver2D(){};
  std::string getName() const { return name; };
  void setGuide( const WaveGuideFDSimulation &guide );
  virtual void build() = 0; // Build the matrices and vectors
protected:
  std::string name;
  const WaveGuideFDSimulation *guide;
};
#endif
