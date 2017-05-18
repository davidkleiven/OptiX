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

  if ( argc != 2 )
  {
    meep::master_printf("No geometry filename specified!");
    return 1;
  }
  Json::Value root;
  Json::Reader reader;
  ifstream infile;
  infile.open(argv[1]);
  if ( !infile.good() )
  {
    clog << "Could not open input file " << argv[1] << endl;
    return 1;
  }
  reader.parse( infile, root );
  infile.close();
  int incidentStokesVector[4]={1,1,0,0};

  srand( time(0) );
  try
  {
    SteadyCoccolithSim sim;
    sim.usePeriodicBoundaryConditions = true;
    sim.prefix = "steadyCoccolith";
    CaCO3Cocco material( 2.19 );

    if ( root.isMember("threshold") )
    {
      material.setThreshold( root["threshold"].asUInt() );
    }
    material.loadRaw( root["voxels"].asString().c_str() );

    if ( root.isMember("incStokesVector") )
    {
      for ( unsigned int i=0;i<4;i++ )
      {
        incidentStokesVector[i] = root["incStokesVector"][i].asInt();
      }
    }
    sim.setIncStokesVector(incidentStokesVector);
    sim.setMaterial( material );
    sim.setMainPropagationDirection( MainPropDirection_t::X );
    sim.setSourceSide( SourcePosition_t::BOTTOM );
    sim.initSource( 0.048 );
    sim.setPMLInWavelengths( 2.0 );

    sim.init();
    sim.BiCStabL = 2; // 2 is default, converges faster with large L, but the memory usage increases
    sim.maxiters = 100000; // 10000 is default
    //sim.run();
    sim.stepForTime( 800.0 );

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
