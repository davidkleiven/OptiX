#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "waveGuideRadiusCurvature.hpp"
#include "cladding.hpp"
#include "numerovSolver.hpp"

class MinimalWaveGuide: public WaveGuide1DSimulation
{
public:
  MinimalWaveGuide():WaveGuide1DSimulation("MinimalGuide"){};
  double potential(double x) const override{ return x;};
};

class NumerovTest: public Numerov
{
public:
  NumerovTest(){eigenvalue=2.0;};
  double alpha_nTest( double x){ return alpha_n(x); };
  double alpha_np1Test( double x ){ return alpha_np1(x); };
  double alpha_nm1Test( double x ){ return alpha_nm1(x); };
  double effPot( double x ){ return effectivePotential(x); };
};

BOOST_AUTO_TEST_SUITE( numerov )
BOOST_AUTO_TEST_CASE( coefficients )
{
  MinimalWaveGuide guide;
  guide.setWavenumber(4.0);

  NumerovTest test;
  test.setStepsize( 0.1 );
  test.setGuide( guide );

  double x = 3.0;
  double alpha = test.alpha_nTest( x );
  double alpha_exp = 241.0/240.0;

  double effectivePot = test.effPot(x);
  double expPot = -(3.0 -2.0);
  BOOST_CHECK_CLOSE( expPot, effectivePot, 0.1f );

  BOOST_CHECK_CLOSE( alpha_exp, alpha, 0.1f );

  alpha = test.alpha_np1Test( x );
  alpha_exp = 1199.0/1200.0;
  BOOST_CHECK_CLOSE( alpha_exp, alpha, 0.1f );

  alpha = test.alpha_nm1Test( x );
  BOOST_CHECK_CLOSE( alpha, alpha_exp, 0.1f );
}
BOOST_AUTO_TEST_SUITE_END()
