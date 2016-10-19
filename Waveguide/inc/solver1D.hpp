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
    void setLimits( double x0, double x2 );
    double getEigenvalue( unsigned int i ) const { return eigenvalues[i]; };
    double getXmin() const { return x1; };
    double getXmax() const { return x2; };
    void addEigenmode( double eigenvalue, double eigvec[], unsigned int size ); // Intended for use when loading
    void loadingFinished(); // Call this after loading is finished
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
    unsigned int nModes{0};
    std::vector<double> eigenvalues;
    arma::mat *solution{NULL};
    const WaveGuide1DSimulation* waveguide{NULL};
    unsigned int maxNumberOfModes{100};

    unsigned int closestIndx( double x ) const;
};
#endif
