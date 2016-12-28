#ifndef SOLVER_H
#define SOLVER_H
#include <string>
#include <vector>
#include <json/writer.h>
#include <armadillo>

class WaveGuide1DSimulation;

/** Base class for the 1D solvers */
class Solver1D
{
  public:
    Solver1D( const char* name ): name(name), solution(new arma::mat()){};
    virtual ~Solver1D();
    virtual void solve() = 0;

    /** Set upper and lower boundaries */
    void setLimits( double x0, double x2 );

    /** Get the i-th eigenvalue (real part) */
    double getEigenvalue( unsigned int i ) const { return eigenvalues[i]; };

    /** Get the imaginary part if the i-th eigenvalue */
    double getEigenvalueImag( unsigned int i ) const { return eigenvaluesImag[i]; };

    /** Get lower boundary */
    double getXmin() const { return x1; };

    /** Get upper boundary */
    double getXmax() const { return x2; };

    /** Add eigenmode. Relevant if loading an old solution */
    void addEigenmode( double eigenvalue, double eigvec[], unsigned int size ); // Intended for use when loading

    /** This function has to be called after an old solution has been loaded */
    void loadingFinished(); // Call this after loading is finished

    /** Get the number of eigenmodes */
    unsigned int getNmodes() const { return nModes; };

    /** Get the number of discretization points */
    unsigned int getEigenVectorSize() const { return solution->n_rows; };

    /** Get the name of the solver */
    std::string getName() const { return name; };

    /** Set waveguide to operate on */
    void setGuide( const WaveGuide1DSimulation &guide );

    /** Get the solution for all eigenmodes */
    const arma::mat& getSolution() const { return *solution; };

    /** Get the solution for all eigenmodes */
    arma::mat& getSolution() { return *solution; };

    /** Get the imaginary part of the solution for all eigenmodes */
    arma::mat& getSolutionImag() { return *solutionImag; };

    /** Get the solution at position x for a given eigenmode */
    double getSolution( double x, unsigned int eigenmode ) const;

    /** Set the number of eigenmodes to store */
    void setNumberOfModesToStore( unsigned int modes ) { nModes=modes;};

    /** Print the solution for the given eigenmode. Only for debugging */
    void printMode( unsigned int mode ) const; // For debugging

    /** Make sure the eigenvectors are normalized */
    double normEigenvec( unsigned int mode ) const;

    /** Fill JSON object with values specific to the class */
    virtual void fillJsonObj( Json::Value &obj ) const;
  protected:
    double x1, x2; // x1: lower limit of domain, x2: upper limit of domain
    std::string name;
    unsigned int nModes{0};
    std::vector<double> eigenvalues;
    std::vector<double> eigenvaluesImag;
    arma::mat *solution{NULL};
    arma::mat *solutionImag{NULL};
    const WaveGuide1DSimulation* waveguide{NULL};
    unsigned int maxNumberOfModes{100};

    /** Get the array index closest to x */
    unsigned int closestIndx( double x ) const;
};
#endif
