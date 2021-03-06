#include "bandStructureSimulation2D.hpp"
#include <iostream>
#include <cassert>
#include <ctime>
#include <iomanip>
#include <cstdlib>
#include <harminv.h>

using namespace std;

BandStructure2DSimulation::~BandStructure2DSimulation()
{
  delete field; field = NULL;
  delete struc; struc = NULL;
  delete srcTime; srcTime = NULL;
  delete ldos; ldos = NULL;
  delete modes; modes= NULL;
}

void BandStructure2DSimulation::printInfo() const
{
  if ( material == NULL )
  {
    throw( runtime_error("No material is set!") );
  }

  if ( meep::am_master() )
  {
    double vxsize = material->info.voxelsize;
    cout << "================= SIMULATION INFORMATION =======================\n";
    cout << "Lx=" << (gdvol.xmax()-gdvol.xmin())*vxsize << " nm. Ly=" << (gdvol.ymax()-gdvol.ymin())*vxsize << " nm.\n";
    cout << "Lambda0=" << vxsize/freq << " nm. lambdaMin=" << vxsize/(freq+freqWidth) << " nm. lambdaMax=" << vxsize/(freq-freqWidth) << " nm\n";
    cout << "Bloch vector: " << bloch.x() << " " << bloch.y() << endl;
    cout << "Numerical resolution: " << vxsize/resolution << " nm\n";
  }
}

void BandStructure2DSimulation::setMaterial( VoxelMaterial2D &newmaterial )
{
  material = &newmaterial;
  gdvol = meep::vol2d( material->Nx(), material->Ny(), resolution );
  material->domain.xmin = gdvol.xmin();
  material->domain.xmax = gdvol.xmax();
  material->domain.ymin = gdvol.ymin();
  material->domain.ymax = gdvol.ymax();
}

void BandStructure2DSimulation::init()
{
  if ( material == NULL )
  {
    throw (runtime_error("No material is set!") );
  }

  printInfo();

  delete struc;
  struc = new meep::structure( gdvol, *material );

  // Update the dispersion
  if ( sellmeier != NULL )
  {
    const double PI = acos(-1.0);
    for ( unsigned int i=0;i<sellmeier->nLorentzians();i++ )
    {
      double sigma, omega0;
      sellmeier->getMEEPLorentzian( material->info.voxelsize*1E-3, i, sigma, omega0 );
      Voxel2DSusceptibility matFunc( sigma, 0.0 );
      struc->add_susceptibility( matFunc, meep::E_stuff, meep::lorentzian_susceptibility(omega0/(2.0*PI), 0.0) );
    }
  }

  delete field;
  field = new meep::fields( struc );
  field->use_bloch( bloch );

  delete srcTime;
  if ( useSingleFrequency )
  {
    srcTime = new meep::continuous_src_time( freq );
  }
  else
  {
    srcTime = new meep::gaussian_src_time( freq-freqWidth/2.0, freq+freqWidth/2.0 );
  }

  // Put the source at more or less a random position
  srcPos = meep::vec( (gdvol.xmax()-gdvol.xmin())*0.3424+gdvol.xmin(), (gdvol.ymax()-gdvol.ymin())*0.64234+gdvol.ymin() );
  unsigned int numberOfTrials = 0;
  unsigned int maxNumberOfTrials = 1000;
  while (( material->isInside(srcPos) ) && (numberOfTrials < maxNumberOfTrials ))
  {
    double randomX = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
    double randomY = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
    srcPos = meep::vec( (gdvol.xmax()-gdvol.xmin())*randomX+gdvol.xmin(), (gdvol.ymax()-gdvol.ymin())*randomY+gdvol.ymin() );
    ++numberOfTrials;
  }

  if ( numberOfTrials == maxNumberOfTrials )
  {
    throw( runtime_error("Did not manage to put the source in a place where it is air!") );
  }
  field->add_point_source( meep::Ez, *srcTime, srcPos );

  // Add DFT point to compute LDOS
  delete ldos;
  ldos = new meep::dft_ldos( freq-freqWidth/2.0, freq+freqWidth/2.0, nfreq );
}

void BandStructure2DSimulation::run()
{
  init();
  double endTime = 10.0*material->Nx();
  if ( !useSingleFrequency )
  {
    endTime = field->last_source_time()+50.0*material->Nx();
  }

  meep::master_printf("End time: %.2f\n", endTime);

  while( field->time() < endTime )
  {
    field->step();
    ldos->update(*field);
    Ez.push_back( field->get_field( meep::Ez, srcPos) );
  }
}

void BandStructure2DSimulation::save()
{
  string fname(prefix);
  if ( addTimestampToFilename )
  {
    setTimeStamp();
    fname += timestamp;
  }

  if ( addUIDToFilename )
  {
    stringstream uidStr;
    uidStr << "_" << uid;
    fname += uidStr.str();
  }
  if ( outfile == NULL )
  {
    outfile = field->open_h5file( fname.c_str() );
  }
  meep::master_printf("Saving parameters...\n");
  saveParams();

  if ( ldos != NULL )
  {
    meep::master_printf("Saving LDOS...\n");
    int size = ldos->Nomega;
    outfile->write("ldos", 1, &size, ldos->ldos(), false );
    outfile->prevent_deadlock();

    size = Ez.size();
    vector<double> EzReal( Ez.size());
    for ( unsigned int i=0;i<Ez.size();i++ )
    {
      EzReal[i] = Ez[i].real();
    }
    meep::master_printf("Saving Ez...\n");
    outfile->write("Ez", 1, &size, &EzReal[0], false);
    outfile->prevent_deadlock();
  }

  if ( modes != NULL )
  {
    meep::master_printf("Saving the HARMINV results...\n");
    outfile->write("freqRe", 1, &modes->numModesFound, modes->freqRe.memptr(), false );
    outfile->prevent_deadlock();

    outfile->write("freqIm", 1, &modes->numModesFound, modes->freqIm.memptr(), false );
    outfile->prevent_deadlock();

    outfile->write("freqErr", 1, &modes->numModesFound, modes->freqErr.memptr(), false );
    outfile->prevent_deadlock();

    // Extract amplitude
    arma::vec amp = arma::abs(modes->amplitude);
    outfile->write("amplitude", 1, &modes->numModesFound, amp.memptr(), false );
    outfile->prevent_deadlock();

    // Extract the phase
    arma::vec phase = arma::arg( modes->amplitude );
    outfile->write("phase", 1, &modes->numModesFound, phase.memptr(), false );
    outfile->prevent_deadlock();
  }
  meep::master_printf("Saving epsilon...\n");
  field->output_hdf5( meep::Dielectric, gdvol.surroundings(), outfile, false, true );

  if ( useSingleFrequency )
  {
    meep::master_printf("Saving spatial Ez...\n");
    field->output_hdf5( meep::Ez, gdvol.surroundings(), outfile, false, true );
  }
  meep::master_printf("All data is saved...\n");
}

void BandStructure2DSimulation::saveParams()
{
  assert( outfile != NULL );

  double PI = acos(-1.0);
  int size = 3;
  double nfreq = ldos->Nomega;
  double freq[size] = {ldos->omega_min/(2.0*PI), ldos->domega/(2.0*PI), nfreq};
  outfile->write("freq", 1, &size, freq, false );
  outfile->prevent_deadlock();

  size = 4;
  double domain[size] = {gdvol.xmin(),gdvol.xmax(),gdvol.ymin(),gdvol.ymax()};
  outfile->write("domain", 1, &size, domain, false );
  outfile->prevent_deadlock();

  size = 2;
  double sourcePosition[size] = {srcPos.x(),srcPos.y()};
  outfile->write("sourcePos", 1, &size, sourcePosition, false );
  outfile->prevent_deadlock();

  double blochVector[size] = {bloch.x(),bloch.y()};
  outfile->write("bloch", 1, &size, blochVector, false );
  outfile->prevent_deadlock();

  size = 1;
  double dt = field->dt;
  outfile->write("dt", 1, &size, &dt, false );
  outfile->prevent_deadlock();
}

void BandStructure2DSimulation::setTimeStamp()
{
  time_t t = time(NULL);
  auto tm = *localtime(&t);
  stringstream ss;
  // Skip seconds as different processes reach this point slightly different, but a minute change in minutes is more unlikely
  ss << put_time( &tm, "%Y%m%d_%H%M");
  timestamp = ss.str();
}

void BandStructure2DSimulation::findModes()
{
  if ( ldos == NULL )
  {
    throw( runtime_error("Find modes require LDOS to extract the max and min frequency!"));
  }
  else if ( field == NULL )
  {
    throw( runtime_error("Find modes require that field is initialized to extract the timestep!"));
  }

  double PI = acos(-1.0);
  modes = new HarminvRes( nHarminvFreq );
  double freqMin = ldos->omega_min/(2.0*PI);
  double freqMax = (ldos->omega_min + ldos->domega*ldos->Nomega)/(2.0*PI);
  modes->numModesFound = meep::do_harminv( &Ez[0], Ez.size(), field->dt, freqMin, freqMax, nHarminvFreq, \
  modes->amplitude.memptr(), modes->freqRe.memptr(), modes->freqIm.memptr(), modes->freqErr.memptr() );
  meep::master_printf("Found %d modes...\n", modes->numModesFound );
}

// ============================ HARMINVRES =====================================
HarminvRes::HarminvRes( unsigned int nfreq )
{
  amplitude.set_size(nfreq);
  freqRe.set_size(nfreq);
  freqIm.set_size(nfreq);
  freqErr.set_size(nfreq);
}
