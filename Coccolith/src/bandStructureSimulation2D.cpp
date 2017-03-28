#include "bandStructureSimulation2D.hpp"
#include <iostream>
#include <cassert>

using namespace std;

BandStructure2DSimulation::~BandStructure2DSimulation()
{
  delete struc; struc = NULL;
  delete field; field = NULL;
  delete srcTime; srcTime = NULL;
  delete ldos; ldos = NULL;
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
    cout << "Lambda0=" << vxsize/freq << " nm. lambdaMin=" << vxsize/(freq+freqWidth) << " nm. lambdaMax=" << vxsize/(freqWidth) << " nm\n";
    cout << "Bloch vector: " << bloch.x() << " " << bloch.y() << endl;
    cout << "Numerical resolution: " << vxsize/resolution << " nm\n";
  }
}

void BandStructure2DSimulation::setMaterial( VoxelMaterial2D &newmaterial )
{
  material = &newmaterial;
  gdvol = meep::vol2d( material->Nx(), material->Ny(), resolution );
}

void BandStructure2DSimulation::init()
{
  if ( material == NULL )
  {
    throw (runtime_error("No material is set!") );
  }

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
  field->add_point_source( meep::Ez, *srcTime, srcPos );

  // Add DFT point to compute LDOS
  delete ldos;
  ldos = new meep::dft_ldos( freq-freqWidth/2.0, freq+freqWidth/2.0, nfreq );
}

void BandStructure2DSimulation::run()
{
  double endTime = 10.0*material->Nx();
  if ( !useSingleFrequency )
  {
    endTime = field->last_source_time()+3.0*material->Nx();
  }

  while( field->time() < endTime )
  {
    field->step();
    ldos->update(*field);
    Ez.push_back( field->get_field( meep::Ez, srcPos).real() );
  }
}

void BandStructure2DSimulation::save( const char* fname )
{
  if ( outfile == NULL )
  {
    outfile = field->open_h5file( fname );
  }
  saveParams();

  if ( ldos != NULL )
  {
    int size = ldos->Nomega;
    outfile->write("ldos", 1, &size, ldos->ldos(), false );
    outfile->prevent_deadlock();

    size = Ez.size();
    outfile->write("Ez", 1, &size, &Ez[0], false );
    outfile->prevent_deadlock();
  }
  field->output_hdf5( meep::Dielectric, gdvol.surroundings(), outfile, false, true );

  if ( useSingleFrequency )
  {
    field->output_hdf5( meep::Ez, gdvol.surroundings(), outfile, false, true );
  }
}

void BandStructure2DSimulation::saveParams()
{
  assert( outfile != NULL );

  int size = 3;
  double nfreq = ldos->Nomega;
  double freq[size] = {ldos->omega_min, ldos->domega, nfreq};
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
}
