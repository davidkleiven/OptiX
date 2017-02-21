#include "voxelMaterial.hpp"
#include "coccolithSim.hpp"
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include<ctime>

using namespace std;

int main( int argc, char** argv )
{
  srand( time(0) );
  try
  {
    //VoxelMaterial material;
    //material.loadRaw( "data/cocco8cv4_216_182_249_253.raw" );
    //material.slideThroughVoxelArray();
    //material.showProjections();

    CoccolithSimulation sim;
    sim.loadVoxels( "data/cocco8cv4_216_182_249_253.raw" );
    sim.setMainPropagationDirection( MainPropDirection_t::X );
    sim.setSourceSide( SourcePosition_t::BOTTOM );
    sim.initSource( 0.25, 0.25);
    sim.setPMLInWavelengths( 2.0 );
    sim.setPlotUpdateFreq( 5 );
    sim.setEndTime( 10.0 );

    sim.runWithoutScatterer();
    sim.init();

    sim.run();
    sim.exportResults();

    sim.runWithScatterer();
    sim.reset();
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
