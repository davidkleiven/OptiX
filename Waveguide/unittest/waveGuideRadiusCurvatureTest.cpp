#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "waveGuideRadiusCurvature.hpp"
#include "cladding.hpp"

BOOST_AUTO_TEST_SUITE( waveGuidRadiusCurv )
BOOST_AUTO_TEST_CASE( potentialTest )
{
  WaveGuideLargeCurvature guide;
  guide.setWidth( 100.0 );
  guide.setWaveLength( 0.1 );
  guide.setRadiusOfCurvature( 40.0E6 );
  Cladding cladding;
  cladding.setElectronDensity( 4066.5 );
  cladding.setRefractiveIndex( 4.14E-5, 3.45E-6 );
  guide.setCladding(cladding);

  double pot = guide.potential(-30.0);
  //double expectedPotential = 0.0059217;
  double expectedPotential = 0.005921;
  BOOST_CHECK_CLOSE( pot, expectedPotential, 0.1f );

  pot = guide.potential(-150.0);
  expectedPotential = 0.3565;
  BOOST_CHECK_CLOSE( pot, expectedPotential, 0.1f );

  pot = guide.potential(150.0);
  expectedPotential = 0.2973;
  BOOST_CHECK_CLOSE( pot, expectedPotential, 0.1f );

}
BOOST_AUTO_TEST_SUITE_END()
