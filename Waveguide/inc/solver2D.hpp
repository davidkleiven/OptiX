#ifndef SOLVER2D_H
#define SOLVER2D_H
#include <string>
#include <complex>

class WaveGuideFDSimulation;

class Solver2D
{
public:
  Solver2D( const char* name ):name(name){};
  virtual ~Solver2D(){};
  std::string getName() const { return name; };
  void setGuide( const WaveGuideFDSimulation &guide );
  virtual void build() = 0; // Build the matrices and vectors
  virtual void solve() = 0;
protected:
  std::string name;
  const WaveGuideFDSimulation *guide;
  std::complex<double> **solution{NULL};

  const std::complex<double>* getSolution( unsigned int iz ) const { return solution[iz]; };
  std::complex<double>* getSolution( unsigned int iz ) { return solution[iz]; };
};
#endif
