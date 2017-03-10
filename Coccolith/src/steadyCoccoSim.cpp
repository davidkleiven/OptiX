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
  centerFrequency = freq;
  contSource = new meep::continuous_src_time( freq );
  sourceTime = contSource;
}

void SteadyCoccolithSim::run()
{
  if ( !isInitialized )
  {
    throw( runtime_error("Simulation is not initialized!") );
  }

  bool success = field->solve_cw( tolerance, maxiters, BiCStabL );
  //bool success = field->solve_cw();


  if ( !success )
  {
    clog << "Convergence failed!\n";
  }
}

void SteadyCoccolithSim::stepForTime( double endTime )
{
  if ( !isInitialized )
  {
    throw( runtime_error("Simulation is not initialized!") );
  }

  while ( field->time() < endTime )
  {
    field->step();
  }
}

void SteadyCoccolithSim::exportResults()
{
  setUID();

  stringstream ss;
  ss << "data/voxelMaterialSimSteady_" << uid; // File extension added automatically
  if ( file == NULL )
  {
    file = field->open_h5file( ss.str().c_str() );
    saveGeometry();
  }
  field->output_hdf5( meep::Dielectric, gdvol.surroundings(), file, false, true );
  file->prevent_deadlock();
  field->output_hdf5( meep::EnergyDensity, gdvol.surroundings(), file, false, true );
  file->prevent_deadlock();
  field->output_hdf5( fieldComp, gdvol.surroundings(), file, false, true );
  file->prevent_deadlock();
  clog << "Results written to " << ss.str() << endl;
}

void SteadyCoccolithSim::init()
{
  if ( material == NULL )
  {
    throw( runtime_error("No material loaded! Call loadVoxels!") );
  }
  domainInfo();
  initializeGeometry();

  material->setDomainSize( gdvol, getPMLThickness() );
  addSourceVolume();
  addStructure();
  addFields();
  addSource();
  isInitialized = true;
}
