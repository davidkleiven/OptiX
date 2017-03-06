#include "voxelMaterial.hpp"
#include "coccolithSim.hpp"
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main( int argc, char** argv )
{
  meep::initialize mpi( argc, argv );
  srand( time(0) );
  try
  {
    CoccolithSimulation sim;
    CaCO3Cocco material;
    material.loadRaw( "data/cocco8cv4Rotated_216_182_249_253.raw" );
    sim.setMaterial( material );
    sim.setMainPropagationDirection( MainPropDirection_t::X );
    sim.setSourceSide( SourcePosition_t::BOTTOM );
    sim.setNfreqFT( 200 );
    sim.initSource( 0.045, 0.03 );
    sim.setPMLInWavelengths( 2.0 );
    sim.setPlotUpdateFreq( 30 );
    sim.disableRealTimeVisualization();
    //sim.setEndTime( 5.0);

    sim.runWithoutScatterer();
    sim.init();

    sim.run();
    sim.exportResults();

    cout << "===============================================================\n";
    cout << "================= RUNNING WITH SCATTERER ======================\n";
    cout << "===============================================================\n";
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

  clog << "Everything has gone out of scope by now. Why is it segfaulting here?\n";
  return 0;
}
