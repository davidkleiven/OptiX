#include "voxelMaterial.hpp"
#include "coccolithSim.hpp"
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <json/reader.h>
#include <fstream>

using namespace std;

int main( int argc, char** argv )
{
  meep::initialize mpi( argc, argv );
  
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

  // Visible light: centerFreq = 0.045, fwidth 0.03 looks good
  double centerFreq = root["centerFreq"].asDouble();
  double freqwidth = root["fwidth"].asDouble();
  bool useDispersive = root["useDispersive"].asBool();
  try
  {
    CoccolithSimulation sim;
    CaCO3Cocco material( 2.72 );
    DispersiveVoxel materialDisp;
    sim.resolution = root["resolution"].asDouble();
    if ( useDispersive )
    {
      materialDisp.loadRaw( "data/cocco8cv4Rotated_216_182_249_253.raw" );
      materialDisp.load( "Materials/CaCO3.json" );
      sim.setMaterial( materialDisp );
    }
    else
    {
      material.loadRaw( "data/cocco8cv4Rotated_216_182_249_253.raw" );
      sim.setMaterial( material );
    }
    sim.setMainPropagationDirection( MainPropDirection_t::X );
    sim.setSourceSide( SourcePosition_t::BOTTOM );
    sim.setNfreqFT( 200 );
    sim.initSource( centerFreq, freqwidth );
    sim.setPMLInWavelengths( 2.0 );
    sim.setPlotUpdateFreq( 30 );
    sim.disableRealTimeVisualization();
    //sim.setEndTime( 5.0);

    sim.runWithoutScatterer();
    sim.init();

    // If the material is dispersive the structure needs to be updated
    if (( useDispersive ) && (meep::my_rank() == 0 ))
    {
      materialDisp.updateStructure( sim.getStructure() );
      clog << "Epsilon at center frequency " << materialDisp.getEpsilon( materialDisp.getVoxelSize()*1E-3, 2.0*3.14159*centerFreq ) << endl;
      clog << "Epsilon at lower frequency " << materialDisp.getEpsilon( materialDisp.getVoxelSize()*1E-3, 2.0*3.14159*(centerFreq-freqwidth)) << endl;
      clog << "Epsilon at upper frequency " <<  materialDisp.getEpsilon( materialDisp.getVoxelSize()*1E-3, 2.0*3.14159*(centerFreq+freqwidth)) << endl;
    }

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
