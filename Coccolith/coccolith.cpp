#include "voxelMaterial.hpp"
#include "coccolithSim.hpp"
#include <stdexcept>
#include <iostream>

using namespace std;

int main( int argc, char** argv )
{
  VoxelMaterial material;
  try
  {
    material.loadRaw( "data/cocco8cv4_216_182_249_253.raw" );
    //material.slideThroughVoxelArray();
    material.showProjections();

    CoccolithSimulation sim;
    sim.loadVoxels( "data/cocco8cv4_216_182_249_253.raw" );
    sim.setMainPropagationDirection( MainPropDirection_t::Y );
    sim.initSource( 2.0, 2.0 );
    sim.setPMLInWavelengths( 0.1 );
    sim.init();
    sim.run();
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
