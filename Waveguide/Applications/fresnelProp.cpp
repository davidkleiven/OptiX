#include "fresnelPropagator.hpp"
#include <string>

class StepSource
{
public:
  double operator()( double x ) const
  {
    if ( ( x > xmin) && ( x < xmax) )
    {
      return 1.0;
    }
    return 0.0;
  }
  double xmin, xmax;
};

using namespace std;

int main( int argc, char** argv )
{
  StepSource source;
  source.xmin = -5.0;
  source.xmax = 5.0;

  double wav = 0.1569;
  FresnelPropagator propagator;
  try
  {
    propagator.setWavelength(wav);
    propagator.setTransverseDiscretization( -100.0, 100.0, 1024 );
    propagator.setInitialConditions( source );
    propagator.setStepsize( 1E3 );
    propagator.propagate( 50 );
    string fname("data/fresnelPropTest.h5");
    propagator.save( fname );
  }
  catch ( exception &exc )
  {
    cout << exc.what() << endl;
    return 1;
  }
  catch (...)
  {
    cout << "Unexpected exception!\n";
    return 1;
  }
  return 0;
}
