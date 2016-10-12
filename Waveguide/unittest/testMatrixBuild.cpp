#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "waveGuideFDSimulation.hpp"
#include "crankNicholson.hpp"
#include <stdexcept>
#include <vector>
#include <armadillo>
#define STEP 0.5

using namespace std;

typedef complex<double> cdouble;

class MinimalWG: public WaveGuideFDSimulation
{
public:
  MinimalWG(){};
  void getXrayMatProp( double x, double z, double &delta, double &beta ) const override final{ delta=_delta; beta=_beta; };
  void setBoundaryConditions() override final;
  cdouble transverseBC( double z ) const override final;
  const Solver2D& getSolver() const { return *solver; };
  double _delta;
  double _beta;
};

cdouble MinimalWG::transverseBC( double z ) const
{
  cdouble im(0.0,1.0);
  return exp(-wavenumber*_beta*z)*exp(-im*_delta*wavenumber*z);
}

void MinimalWG::setBoundaryConditions()
{
  vector<cdouble> values(nodeNumberTransverse(), 1.0);
  solver->setLeftBC(&values[0]);
};

MinimalWG* initWG()
{
  cdouble imgu(0.0,1.0);
  MinimalWG *wg = new MinimalWG();
  wg->setTransverseDiscretization(0.0, 10.0, STEP);
  wg->setLongitudinalDiscretization(0.0,10.0, STEP);
  wg->setWavenumber( 0.1 );
  return wg;
}
BOOST_AUTO_TEST_SUITE( testSolver )
BOOST_AUTO_TEST_CASE( testBCSettings )
{

  MinimalWG *wg = initWG();
  CrankNicholson solver;
  wg->setSolver(solver);
  wg->setBoundaryConditions();

  unsigned int Nx = wg->nodeNumberTransverse();
  unsigned int Nz = wg->nodeNumberLongitudinal();
  const arma::cx_mat bc = wg->getSolver().getSolution(0);
  unsigned int bcUnequal = 0;
  for ( unsigned int ix=0;ix<Nx;ix++ )
  {
    //cout << bc(ix,0) << endl;
    if ( abs(bc(ix,0).real() - 1.0) > 1E-6 )
    {
      bcUnequal++;
    }
  }
  delete wg;
  BOOST_CHECK_EQUAL( bcUnequal, 0 );
}

BOOST_AUTO_TEST_CASE( testSolution )
{
  MinimalWG *wg = initWG();
  CrankNicholson solver;
  wg->setSolver(solver);
  wg->setBoundaryConditions();

  unsigned int Nx = wg->nodeNumberTransverse();
  unsigned int Nz = wg->nodeNumberLongitudinal();
  const arma::cx_mat bc = wg->getSolver().getSolution(0);

  wg->solve();

  const arma::cx_mat sol = wg->getSolver().getSolution(0);

  double expectedSolution=1.0;

  unsigned int numberOfUnequalElements = 0;
  for ( unsigned int ix=0;ix<Nx; ix++ )
  {
    for ( unsigned int iz=0;iz<Nz; iz++ )
    {
      if ( abs( sol(ix,iz).real() - expectedSolution) > 1E-4 )
      {
        numberOfUnequalElements++;
      }
    }
  }
  delete wg;
  BOOST_CHECK_EQUAL( numberOfUnequalElements, 0 );
}

BOOST_AUTO_TEST_CASE( testWidthAbs )
{
  MinimalWG *wg = initWG();
  CrankNicholson solver;
  wg->setSolver(solver);
  wg->setBoundaryConditions();
  wg->_beta = 0.1;

  unsigned int Nx = wg->nodeNumberTransverse();
  unsigned int Nz = wg->nodeNumberLongitudinal();
  const arma::cx_mat bc = wg->getSolver().getSolution(0);

  wg->solve();
  unsigned int numberOfUnequalElements = 0;
  const arma::cx_mat sol = wg->getSolver().getSolution(0);
  for ( unsigned int ix=0;ix<Nx; ix++ )
  {
    for ( unsigned int iz=0;iz<Nz; iz++ )
    {
      double z = iz*STEP;
      double expected = exp(-wg->_beta*z*0.1);
      if ( abs( sol(ix,iz).real() - expected) > 1E-2 )
      {
        numberOfUnequalElements++;
      }
    }
  }
  delete wg;
  BOOST_CHECK_EQUAL( numberOfUnequalElements, 0 );
}
BOOST_AUTO_TEST_SUITE_END()
