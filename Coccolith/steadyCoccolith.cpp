#include "voxelMaterial.hpp"
#include "steadyCoccoSim.hpp"
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main( int argc, char** argv )
{
  srand( time(0) );
  try
  {
    SteadyCoccolithSim sim;
    sim.loadVoxels( "data/cocco8cv4Rotated_216_182_249_253.raw" );
    sim.setMainPropagationDirection( MainPropDirection_t::X );
    sim.setSourceSide( SourcePosition_t::BOTTOM );
    sim.initSource( 0.036 );
    sim.setPMLInWavelengths( 2.0 );
    sim.init();
    sim.run();
    sim.exportResults();
  }
  catch( exception &exc )
  {
    cout << exc.what() << endl;
  }
  catch( ... )
  {
    cout << "Unrecognized exception!\n";
  }

  clog << "Everything has gone out of scope by now. Why is it segfaulting here?\n";
  return 0;
}
