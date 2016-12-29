#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "curvedWaveGuide2D.hpp"

class CWGTest: public CurvedWaveGuideFD
{
public:
  bool isIn( double x, double z ) const { return isInsideGuide(x,z);};
};

BOOST_AUTO_TEST_SUITE( curvedWG )
BOOST_AUTO_TEST_CASE( insideWG )
{
  CWGTest wg;
  wg.setRadiusOfCurvature( 40E6 );
  wg.setWidth( 100.0 );
  bool inside = wg.isIn( 50.0, 0.0 );
  BOOST_CHECK_EQUAL(inside, true);

  inside = wg.isIn(-1.0,0.0);
  BOOST_CHECK_EQUAL(inside, false);

  inside = wg.isIn(0.0, 100.0);
  BOOST_CHECK_EQUAL( inside, true);
}
BOOST_AUTO_TEST_SUITE_END()
