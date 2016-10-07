#ifndef SOLVER_H
#define SOLVER_H
#include <string>
#include <vector>
#include <jsoncpp/json/writer.h>

class WaveGuide1DSimulation;
class Solver1D
{
  public:
    Solver1D( const char* name ): name(name), solution(new std::vector<double>){};
    virtual ~Solver1D();
    virtual void solve() = 0;
    double getEigenvalue() const { return eigenvalue; };
    std::string getName() const { return name; };
    void setGuide( const WaveGuide1DSimulation &guide ){ waveguide = &guide; };
    const std::vector<double>& getSolution() const { return *solution; };
    virtual void fillJsonObj( Json::Value &obj ) const;
  protected:
    std::string name;
    double eigenvalue;
    std::vector<double> *solution;
    const WaveGuide1DSimulation* waveguide{NULL};
};
#endif
