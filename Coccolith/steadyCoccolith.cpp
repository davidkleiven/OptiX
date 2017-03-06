#include "voxelMaterial.hpp"
#include "steadyCoccoSim.hpp"
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <mpi.h>

using namespace std;

int main( int argc, char** argv )
{
  meep::initialize mpi( argc, argv );
  srand( time(0) );
  try
  {
    SteadyCoccolithSim sim;
    CaCO3Cocco material;
    material.loadRaw( "data/cocco8cv4Rotated_216_182_249_253.raw" );
    sim.setMaterial( material );
    sim.setMainPropagationDirection( MainPropDirection_t::X );
    sim.setSourceSide( SourcePosition_t::BOTTOM );
    sim.initSource( 0.036 );
    sim.setPMLInWavelengths( 2.0 );
    sim.init();
    sim.BiCStabL = 2; // 2 is default, converges faster with large L, but the memory usage increases
    sim.maxiters = 100000; // 10000 is default
    //sim.run();
    sim.stepForTime( 100.0 );
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

  clog << "Process " << meep::my_rank() << " finished!\n";
  return 0;
}
