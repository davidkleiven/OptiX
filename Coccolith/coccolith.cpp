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
    sim.setNfreqFT( 50 );
    sim.initSource( 0.25, 0.25 );
    sim.setPMLInWavelengths( 2.0 );
    sim.setPlotUpdateFreq( 30 );
    sim.setEndTime( 5.0);

//    sim.runWithoutScatterer();
    sim.init();

    sim.run();
    /*
    sim.exportResults();

    cout << "===============================================================\n";
    cout << "================= RUNNING WITH SCATTERER ======================\n";
    cout << "===============================================================\n";
    sim.runWithScatterer();
    sim.reset();
    sim.run();
    */sim.exportResults();
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
