#include "steadyCoccoSim.hpp"
#include <stdexcept>
#include <iostream>

using namespace std;

SteadyCoccolithSim::~SteadyCoccolithSim()
{
  if ( contSource != NULL ) delete contSource;
}

void SteadyCoccolithSim::initSource( double freq )
{
  if ( geoIsInitialized )
  {
    throw( runtime_error("Geometry is already initialized! PML layers may behave strangely if the frequencies change!") );
  }

  if ( contSource != NULL ) delete contSource;
  contSource = new meep::continuous_src_time( freq );
  sourceTime = contSource;
}

void SteadyCoccolithSim::run()
{
  if ( !isInitialized )
  {
    throw( runtime_error("Simulation is not initialized!") );
  }

  field->solve_cw();
}

void SteadyCoccolithSim::exportResults()
{
  if ( uid == 0 )
  {
    uid = rand()%UID_MAX;
  }

  stringstream ss;
  ss << "data/voxelMaterialSimSteady_" << uid; // File extension added automatically
  if ( file == NULL )
  {
    file = field->open_h5file( ss.str().c_str() );
    saveGeometry();
  }
  field->output_hdf5( meep::Dielectric, gdvol.surroundings(), file, false, true );
  field->output_hdf5( meep::EnergyDensity, gdvol.surroundings(), file, false, true );
  clog << "Results written to " << ss.str() << endl;
}

void SteadyCoccolithSim::init()
{
  if ( !materialLoaded )
  {
    throw( runtime_error("No material loaded! Call loadVoxels!") );
  }
  domainInfo();
  initializeGeometry();

  material.setDomainSize( gdvol, getPMLThickness() );
  addSourceVolume();
  addStructure();
  addFields();
  addSource();
  isInitialized = true;
}
