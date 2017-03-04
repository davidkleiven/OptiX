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
    sim.BiCStabL = 2; // 2 is default, converges faster with large L, but the memory usage increases
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
  return 0;
}
