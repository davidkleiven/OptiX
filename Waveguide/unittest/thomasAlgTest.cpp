#include <boost/test/unit_test.hpp>
#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN
#include "thomasAlgorithm.hpp"
BOOST_AUTO_TEST_SUITE( ThomasAlgo )

BOOST_AUTO_TEST_CASE( ThomasTest )
{
  ThomasAlgorithm solver;
  double diag[3] = {2.0, 4.0, 5.0};
  double subdiag[2] = {-1.0, 2.0};
  double rhs[3] = {1.0, 2.0, 3.0};
  solver.solve(diag, subdiag, rhs, 3);

  double solution[3]={ 20.0/27.0, 13.0/27.0, 11.0/27.0};
  for ( unsigned int i=0;i<3;i++ )
  {
    BOOST_CHECK_CLOSE(diag[i], solution[i], 0.1f);
  }
}

BOOST_AUTO_TEST_SUITE_END()
