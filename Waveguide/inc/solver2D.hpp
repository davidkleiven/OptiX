#ifndef SOLVER2D_H
#define SOLVER2D_H
#include <string>
#include <complex>
#include <jsoncpp/json/writer.h>
#include <armadillo>

class WaveGuideFDSimulation;

typedef std::complex<double> cdouble;


class Solver2D
{
public:
  enum class Comp_t{REAL,IMAG};
  Solver2D( const char* name ):name(name){};
  virtual ~Solver2D();
  std::string getName() const { return name; };
  void setGuide( const WaveGuideFDSimulation &guide );
  const arma::cx_mat& getSolution( unsigned int iz ) const { return *solution; }; // Depricated. iz is not used.
  const arma::cx_mat& getSolution() const { return *solution; };
  bool importHDF5( const std::string &fname );
  bool importHDF5( const std::string &realpart, const std::string &amplitude );
  void realPart( double *realsolution ) const;
  void imagPart( double *imagsolution ) const;

  // Set boundary condition at z=0
  void setLeftBC( const cdouble values[] );

  // Virtual functions
  virtual void solve() = 0;
  virtual void fillInfo( Json::Value &obj ) const;
protected:
  std::string name;
  const WaveGuideFDSimulation *guide;
  arma::cx_mat *solution{NULL};

  arma::cx_mat& getSolution( unsigned int iz ) { return *solution; };
  void realOrImagPart( double *solution, Comp_t comp ) const;
};
#endif
