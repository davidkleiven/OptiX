#include <meep.hpp>
#include <iostream>
#include <complex>

using namespace std;
typedef complex<double> cdouble;

// Global variables
meep::h5file *file = NULL;
double radius = 10.0;
unsigned int numberOfFreq = 400;
volatile bool overwriteFile = true;
double freq = 0.4, fwidth = 0.2;
double pmlThick = 3.0/freq;
double Nx = 2.0*(radius+pmlThick);
double Ny = 2.0*(radius+pmlThick);
double Nz = 2.0*(radius+pmlThick);

double eps( const meep::vec&r )
{
  if ( pow(r.x()-Nx/2.0,2) + pow(r.y()-Ny/2.0,2) + pow(r.z()-Nz/2.0,2) < radius*radius )
  {
    //return 3.0;
    return -1E12;
  }
  return 1.0;
}

cdouble amplitude( const meep::vec &r )
{
  return 1.0;
}

double epsRef( const meep::vec &r )
{
  return 1.0;
}

void run( double (*eps)(const meep::vec &r) )
{
  double resolution = 7.0; // pixels per distance
  meep::grid_volume v = meep::vol3d(Nx,Ny,Nz, resolution); // 5x10 2d cell
  meep::structure struc(v, eps, meep::pml(pmlThick));
  meep::fields f(&struc);
  meep::gaussian_src_time src(freq, fwidth);

  meep::volume srcVol( meep::vec(pmlThick,0.0,0.0), meep::vec(pmlThick,Ny,Nz) );
  f.add_volume_source( meep::Ey, src, srcVol, amplitude );

  // Add flux planes
  meep::vec refCrn1( Nx/2.0-radius-1.0, Ny/2.0-radius-1.0, Nz/2.0 - radius- 1.0 );
  meep::vec refCrn2( Nx/2.0-radius-1.0, Ny/2.0+radius+1.0, Nz/2.0+radius+1.0);
  meep::volume refPlaneVol(refCrn1, refCrn2);
  meep::vec boxEndCrn( Nx/2.0+radius+1.0, Ny/2.0+radius+1.0, Nz/2.0+radius+1.0);
  meep::volume fluxBoxVol( refCrn1, boxEndCrn );
  meep::dft_flux refPlane = f.add_dft_flux_plane( refPlaneVol, freq-fwidth/2.0, freq+fwidth/2.0, numberOfFreq );
  meep::dft_flux fluxBox = f.add_dft_flux_box( fluxBoxVol, freq-fwidth/2.0, freq+fwidth/2.0, numberOfFreq );

  if (!overwriteFile)
  {
    // This is the second time the funciton is called
    refPlane.load_hdf5(f, "data/refPlaneSphereRef");
    fluxBox.load_hdf5(f, "data/fluxBoxSphereRef");
    refPlane.scale_dfts(-1.0);
    fluxBox.scale_dfts(-1.0);
  }


  double timeToProp = 10.0*Nz;
  meep::master_printf("Source last time %.2f\n", src.last_time() );
  meep::master_printf("Time to propagate %.2f\n", timeToProp );
  meep::master_printf("End time: %.2f\n", src.last_time()+timeToProp);
  while ( f.time() < src.last_time() + timeToProp )
  {
    f.step();
  }

  if ( overwriteFile )
  {
    refPlane.save_hdf5( f, "data/refPlaneSphereRef" );
    fluxBox.save_hdf5( f, "data/fluxBoxSphereRef" );
    file = f.open_h5file("data/sphereCrossSection");

    int length = 3;
    double freqInfo[length];
    freqInfo[0] = refPlane.freq_min;
    freqInfo[1] = refPlane.Nfreq;
    freqInfo[2] = refPlane.dfreq;
    file->write( "freq", 1, &length, freqInfo, false );
    file->prevent_deadlock();

    double *fluxRef = refPlane.flux();
    int Nfreq = refPlane.Nfreq;
    file->write("refPlaneFlux", 1, &Nfreq, fluxRef, false );
    file->prevent_deadlock();

    double *fluxOut = fluxBox.flux();
    file->write("refFluxBox", 1, &Nfreq, fluxOut, false );
    file->prevent_deadlock();

    int nCrn = 6;
    meep::vec minCrn = refPlaneVol.get_min_corner();
    meep::vec maxCrn = refPlaneVol.get_max_corner();
    double crnRefPlane[nCrn] = {minCrn.x(), minCrn.y(), minCrn.z(), maxCrn.x(), maxCrn.y(), maxCrn.z()};
    file->write( "refPlaneCrd", 1, &nCrn, crnRefPlane, false );
    file->prevent_deadlock();
  }
  else
  {
    file = f.open_h5file("data/sphereCrossSection", meep::h5file::READWRITE );
    double *fluxRef = refPlane.flux();
    int Nfreq = refPlane.Nfreq;
    file->write("scatPlaneFlux", 1, &Nfreq, fluxRef, false );
    file->prevent_deadlock();

    double *fluxOut = fluxBox.flux();
    file->write("scatFluxBox", 1, &Nfreq, fluxOut, false );
    file->prevent_deadlock();
    f.output_hdf5( meep::Dielectric, v.surroundings(), file, false, true );
    file->prevent_deadlock();
  }
}
int main(int argc, char **argv) {
  meep::initialize mpi(argc, argv); // do this even for non-MPI Meep
  try
  {
    run( epsRef );
    overwriteFile = false;
    run( eps );
  }
  catch(...)
  {
    clog << "An exception occured...\n!";
  }
  return 0;
}
