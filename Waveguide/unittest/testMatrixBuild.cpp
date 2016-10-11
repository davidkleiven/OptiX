#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "waveGuideFDSimulation.hpp"
#include "crankNicholson.hpp"
#include <stdexcept>

using namespace std;

typedef complex<double> cdouble;

class MinimalWG: public WaveGuideFDSimulation
{
public:
  MinimalWG();
  cdouble getRefractiveIndex( double x, double z) const override final{ return 1.0; };
  void setBoundaryConditions() override final;
};

MinimalWG::MinimalWG()
{
  setTransverseDiscretization(0.0, 2.0, 1.0);
  setLongitudinalDiscretization(0.0, 2.0, 1.0);
  setWavenumber(1.0);
}

void MinimalWG::setBoundaryConditions()
{
  cdouble bc[3] = {0.0,1.0,0.0};
  solver->setLeftBC(bc);
};

class CNTest: public CrankNicholson
{
public:
  CNTest(){};
  void solveSingle( unsigned int iz){ return solveCurrent(iz); };
};

BOOST_AUTO_TEST_CASE( testMatrixSolve )
{
  cdouble imgu(0.0,1.0);
  MinimalWG wg;
  CNTest cn;
  try
  {
    wg.setSolver( cn );
    wg.setBoundaryConditions();
    cn.initValuesFromWaveGuide();
  }
  catch ( exception &exc )
  {
    cout << exc.what() << endl;
    return;
  }
  BOOST_REQUIRE( wg.nodeNumberTransverse() == 3 );
  BOOST_REQUIRE( wg.nodeNumberLongitudinal() == 3 );

  ThomasAlgorithm ta;
  cn.solveSingle( 1 );
  cdouble alpha = imgu/2.0;
  cdouble gamma = alpha*wg.getRefractiveIndex(0.0,0.0);

  cdouble diag[3] = {1.0+alpha-gamma/2.0, 1.0+alpha-gamma/2.0, 1.0+alpha-gamma/2.0};
  cdouble subdiag[2] = {-0.5*alpha, -0.5*alpha};
  cdouble rhs[3] = { 0.5*alpha, -alpha+1.0+0.5*gamma, 0.5*alpha};
  ta.solve( diag, subdiag, rhs, 3);
  const cdouble* sol = static_cast<const CNTest&>(cn).getSolution(1);

  for ( unsigned int i=0;i<3;i++ )
  {
    BOOST_CHECK_CLOSE( sol[i].real(), diag[i].real(), 0.1f );
    BOOST_CHECK_CLOSE( sol[i].imag(), diag[i].imag(), 0.1f );
  }
}
