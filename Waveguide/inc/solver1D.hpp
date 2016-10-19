#ifndef SOLVER_H
#define SOLVER_H
#include <string>
#include <vector>
#include <jsoncpp/json/writer.h>
#include <armadillo>

class WaveGuide1DSimulation;
class Solver1D
{
  public:
    Solver1D( const char* name ): name(name), solution(new arma::mat()){};
    virtual ~Solver1D();
    virtual void solve() = 0;
    double getEigenvalue( unsigned int i ) const { return eigenvalues[i]; };
    unsigned int getNmodes() const { return nModes; };
    unsigned int getEigenVectorSize() const { return solution->n_rows; };
    std::string getName() const { return name; };
    void setGuide( const WaveGuide1DSimulation &guide ){ waveguide = &guide; };
    const arma::mat& getSolution() const { return *solution; };
    void setNumberOfModesToStore( unsigned int modes ) { nModes=modes;};
    virtual void fillJsonObj( Json::Value &obj ) const;
  protected:
    std::string name;
    unsigned int nModes{1};
    std::vector<double> eigenvalues;
    arma::mat *solution{NULL};
    const WaveGuide1DSimulation* waveguide{NULL};
};
#endif
