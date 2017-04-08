#include "voxelMaterial.hpp"
#include "coccolithSim.hpp"
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <json/reader.h>
#include <fstream>
#define CHECK_THAT_ALL_WORKS

using namespace std;

void setupSim( CoccolithSimulation &sim )
{

}

int main( int argc, char** argv )
{
  srand(time(0));
  try
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

    int incidentStokesVector[4]={1,1,0,0};

    // Visible light: centerFreq = 0.045, fwidth 0.03 looks good
    double centerFreq = root["centerFreq"].asDouble();
    double freqwidth = root["fwidth"].asDouble();
    bool useDispersive = root["useDispersive"].asBool();
    CoccolithSimulation *sim = new CoccolithSimulation();
    sim->prefix = root["prefix"].asString();
    SellmeierMaterial sellmeier;
    sim->resolution = root["resolution"].asDouble();
    sellmeier.load( root["material"].asString().c_str() );
    VoxelSusceptibility material( sellmeier.epsInf, 1.0 );
    material.loadRaw( root["voxels"].asString().c_str() );
    sim->setMaterial( material );

    unsigned int nFreq = 150;
    double pmlThick = 2.0;
    sim->setMainPropagationDirection( MainPropDirection_t::X );
    sim->setSourceSide( SourcePosition_t::BOTTOM );
    sim->setNfreqFT( nFreq );
    sim->initSource( centerFreq, freqwidth );
    sim->setPMLInWavelengths( pmlThick );
    sim->disableRealTimeVisualization();

    #ifdef CHECK_THAT_ALL_WORKS
      sim->setEndTime( 10.0);
    #endif

    sim->additionalVaccumLayerPx = 3.0;
    if ( root.isMember("computeAsymmetryFactor") ) sim->computeAsymmetryFactor = root["computeAsymmetryFactor"].asBool();
    if ( root.isMember("computeStokes") ) sim->computeStokesParameters = root["computeStokes"].asBool();
    if ( root.isMember("incStokesVector") )
    {
      for ( unsigned int i=0;i<4;i++ )
      {
        incidentStokesVector[i] = root["incStokesVector"][i].asInt();
      }
    }
    sim->setIncStokesVector(incidentStokesVector);
    sim->runWithoutScatterer();
    sim->init();
    sim->run();
    sim->exportResults();
    string uid = sim->uid;
    delete sim;

    sim = new CoccolithSimulation();
    sim->prefix = root["prefix"].asString();
    sim->resolution = root["resolution"].asDouble();
    sim->uid = uid;
    sim->gaussLegendreOrder = 128;
    sim->numberOfAzimuthalSteps = 100;
    if ( useDispersive )
    {
      sim->setSellmeierMaterial( sellmeier );
    }
    if ( root.isMember("threshold") )
    {
      material.setThreshold( root["threshold"].asUInt() );
    }

    sim->setIncStokesVector(incidentStokesVector);
    sim->setMaterial( material );
    sim->setMainPropagationDirection( MainPropDirection_t::X );
    sim->setSourceSide( SourcePosition_t::BOTTOM );
    sim->setNfreqFT( nFreq );
    sim->initSource( centerFreq, freqwidth );
    sim->setPMLInWavelengths( pmlThick );
    sim->disableRealTimeVisualization();

    #ifdef CHECK_THAT_ALL_WORKS
      sim->setEndTime( 10.0 );
      sim->gaussLegendreOrder = 2;
    #endif

    sim->additionalVaccumLayerPx = 3.0;
    if ( root.isMember("computeAsymmetryFactor") ) sim->computeAsymmetryFactor = root["computeAsymmetryFactor"].asBool();
    if ( root.isMember("computeStokes") ) sim->computeStokesParameters = root["computeStokes"].asBool();
    sim->runWithScatterer();
    sim->init();
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
