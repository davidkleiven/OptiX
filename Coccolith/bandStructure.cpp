#include "bandStructureSimulation2D.hpp"
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <mpi.h>
#include <json/json.h>

using namespace std;

int main( int argc, char** argv )
{
  meep::initialize mpi( argc, argv );
  if ( argc != 2 )
  {
    meep::master_printf("No inputfile specified!");
    return 1;
  }
  string inputfname(argv[1]);
  try
  {
    Json::Value root;
    Json::Reader reader;
    ifstream infile;
    infile.open(inputfname.c_str());
    if ( !infile.good() )
    {
      meep::master_printf("Could not open input file!\n");
      return 1;
    }
    reader.parse( infile, root );
    infile.close();

    SellmeierMaterial material;
    material.load( root["material"].asString().c_str() );

    Voxel2DSusceptibility geom( material.epsInf, 1.0 );
    geom.threshold = 10;
    geom.loadRaw( root["geometry"].asString() );

    BandStructure2DSimulation sim;
    sim.prefix = root["prefix"].asString();
    sim.uid = root["uid"].asUInt();
    sim.addUIDToFilename = root["appendUIDToFileName"].asBool();
    sim.addTimestampToFilename = root["timestampInFilename"].asBool();
    sim.resolution = root["resolution"].asDouble();
    sim.freq = root["freq"].asDouble();
    sim.freqWidth = root["freqWidth"].asDouble();
    sim.nfreq = root["nfreq"].asUInt();
    sim.bloch = meep::vec( root["bloch"][0].asDouble(), root["bloch"][1].asDouble() );
    sim.sellmeier = &material;
    sim.setMaterial( geom );
    sim.run();
    sim.save();
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
