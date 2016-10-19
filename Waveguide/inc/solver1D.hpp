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
    void addEigenvalue( double eigval ) { eigenvalues.push_back(eigval); }; // Intended for use when loading
    unsigned int getNmodes() const { return nModes; };
    unsigned int getEigenVectorSize() const { return solution->n_rows; };
    std::string getName() const { return name; };
    void setGuide( const WaveGuide1DSimulation &guide ){ waveguide = &guide; };
    const arma::mat& getSolution() const { return *solution; };
    arma::mat& getSolution() { return *solution; };
    double getSolution( double x, unsigned int eigenmode ) const;
    void setNumberOfModesToStore( unsigned int modes ) { nModes=modes;};
    virtual void fillJsonObj( Json::Value &obj ) const;
  protected:
    double x1, x2; // x1: lower limit of domain, x2: upper limit of domain
    std::string name;
    unsigned int nModes{1};
    std::vector<double> eigenvalues;
    arma::mat *solution{NULL};
    const WaveGuide1DSimulation* waveguide{NULL};

    unsigned int closestIndx( double x ) const;
};
#endif
