#include "voxelMaterial.hpp"
#include "coccolithSim.hpp"
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <json/reader.h>
#include <fstream>

using namespace std;

void setupSim( CoccolithSimulation &sim )
{

}

int main( int argc, char** argv )
{
  try
  {
    meep::initialize mpi( argc, argv );
    DispersiveVoxel materialDisp;
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
    CoccolithSimulation *sim = new CoccolithSimulation();
    CaCO3Cocco material( 2.72 );
    sim->resolution = root["resolution"].asDouble();
    if ( useDispersive )
    {
      materialDisp.loadRaw( "data/cocco8cv4Rotated_216_182_249_253.raw" );
      materialDisp.load( "Materials/CaCO3.json" );
      sim->setMaterial( materialDisp );
    }
    else
    {
      material.loadRaw( "data/cocco8cv4Rotated_216_182_249_253.raw" );
      sim->setMaterial( material );
    }
    unsigned int nFreq = 200;
    double pmlThick = 2.0;
    sim->setMainPropagationDirection( MainPropDirection_t::X );
    sim->setSourceSide( SourcePosition_t::BOTTOM );
    sim->setNfreqFT( nFreq );
    sim->initSource( centerFreq, freqwidth );
    sim->setPMLInWavelengths( pmlThick );
    sim->disableRealTimeVisualization();
    //sim->setEndTime( 10.0);
    sim->runWithoutScatterer();
    sim->init();
    sim->run();
    sim->exportResults();
    string uid = sim->uid;
    delete sim;

    sim = new CoccolithSimulation();
    sim->resolution = root["resolution"].asDouble();
    sim->uid = uid;
    if ( useDispersive )
    {
      sim->setMaterial( materialDisp );
    }
    else
    {
      sim->setMaterial( material );
    }
    sim->setMainPropagationDirection( MainPropDirection_t::X );
    sim->setSourceSide( SourcePosition_t::BOTTOM );
    sim->setNfreqFT( nFreq );
    sim->initSource( centerFreq, freqwidth );
    sim->setPMLInWavelengths( pmlThick );
    sim->disableRealTimeVisualization();
    //sim->setEndTime( 10.0);
    sim->runWithScatterer();
    sim->init();
    if ( useDispersive )
    {
      materialDisp.updateStructure( sim->getStructure() );
    }
    sim->run();
    sim->exportResults();
    clog << "Process " << meep::my_rank() << " finished\n";
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
