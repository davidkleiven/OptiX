#include <meep.hpp>
#include <iostream>
#include <complex>
#define GAUSSIAN_SOURCE

using namespace std;
typedef complex<double> cdouble;

// Global variables
meep::h5file *file = NULL;
double radius = 3.0;
unsigned int numberOfFreq = 50;
volatile bool overwriteFile = true;
double freq = 0.4, fwidth = 0.3;
double pmlThick = 3.0/freq;
double Nx = 2.0*(1.5*radius+pmlThick);
double Ny = 2.0*(1.5*radius+pmlThick);
double Nz = 2.0*(1.5*radius+pmlThick);

double eps( const meep::vec&r )
{
  if ( pow(r.x()-Nx/2.0,2) + pow(r.y()-Ny/2.0,2) + pow(r.z()-Nz/2.0,2) < radius*radius )
  {
    return 2.0;
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

void run( double (*epsilonFunc)(const meep::vec &r) )
{
  double resolution = 16.0; // pixels per distance
  meep::grid_volume v = meep::vol3d(Nx,Ny,Nz, resolution); // 5x10 2d cell
  meep::structure struc(v, epsilonFunc, meep::pml(pmlThick));
  meep::fields f(&struc);

  #ifdef GAUSSIAN_SOURCE
    meep::gaussian_src_time src(freq, fwidth);
  #else
    meep::continuous_src_time src(freq);
  #endif

  meep::volume srcVol( meep::vec(pmlThick,0.0,0.0), meep::vec(pmlThick,Ny,Nz) );
  meep::component fieldComp = meep::Ez;
  f.add_volume_source( fieldComp, src, srcVol, amplitude );

  #ifdef GAUSSIAN_SOURCE
    // Add flux planes
    //meep::vec refCrn1( Nx/2.0-radius-1.0, Ny/2.0-radius-1.0, Nz/2.0-radius-1.0 );
    //meep::vec refCrn2( Nx/2.0-radius-1.0, Ny/2.0+radius+1.0, Nz/2.0+radius+1.0);
    //meep::volume refPlaneVol(refCrn1, refCrn2);
    //meep::vec boxEndCrn( Nx/2.0+radius+1.0, Ny/2.0+radius+1.0, Nz/2.0+radius+1.0);
    //meep::vec transPlaneCrn( Nx/2.0+radius+1.0, Ny/2.0-radius-1.0, Nz/2.0-radius-1.0 );
    //meep::volume fluxBoxVol( refCrn1, boxEndCrn );
    //meep::volume transPlaneVol( transPlaneCrn, boxEndCrn );
    //meep::dft_flux refPlane = f.add_dft_flux_plane( refPlaneVol, freq-fwidth/2.0, freq+fwidth/2.0, numberOfFreq );
    //meep::dft_flux transPlane = f.add_dft_flux_plane( transPlaneVol, freq-fwidth/2.0, freq+fwidth/2.0, numberOfFreq );
    //meep::dft_flux fluxBox = f.add_dft_flux_box( fluxBoxVol, freq-fwidth/2.0, freq+fwidth/2.0, numberOfFreq );

    meep::vec refCrn1( pmlThick+2.0, pmlThick, pmlThick );
    meep::vec refCrn2( pmlThick+2.0, Ny-pmlThick-1.0, Nz-pmlThick-1.0);
    meep::volume refPlaneVol(refCrn1, refCrn2);
    meep::vec boxEndCrn( Nx-pmlThick-1.0, Ny-pmlThick-1.0, Nz-pmlThick-1.0);
    meep::vec transPlaneCrn( Nx-pmlThick-1.0, pmlThick, pmlThick );
    meep::volume fluxBoxVol( refCrn1, boxEndCrn );
    meep::volume transPlaneVol( transPlaneCrn, boxEndCrn );
    meep::dft_flux refPlane = f.add_dft_flux_plane( refPlaneVol, freq-fwidth/2.0, freq+fwidth/2.0, numberOfFreq );
    meep::dft_flux transPlane = f.add_dft_flux_plane( transPlaneVol, freq-fwidth/2.0, freq+fwidth/2.0, numberOfFreq );
    meep::dft_flux fluxBox = f.add_dft_flux_box( fluxBoxVol, freq-fwidth/2.0, freq+fwidth/2.0, numberOfFreq );

    if (!overwriteFile)
    {
      // This is the second time the function is called
      refPlane.load_hdf5(f, "data/refPlaneSphereRef");
      fluxBox.load_hdf5(f, "data/fluxBoxSphereRef");
      refPlane.scale_dfts(-1.0);
      fluxBox.scale_dfts(-1.0);
    }
  #endif


  double timeToProp = 4.0*Nz;
  #ifdef GAUSSIAN_SOURCE
    double endTime = src.last_time() + timeToProp;
    meep::master_printf("Source last time %.2f\n", src.last_time() );
    meep::master_printf("Time to propagate %.2f\n", timeToProp );
    meep::master_printf("End time: %.2f\n", endTime);
  #else
    double endTime = 50.0*Nz;
  #endif

  while ( f.time() < endTime)
  {
    f.step();
  }

  #ifdef GAUSSIAN_SOURCE
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

      double *fluxTrans = transPlane.flux();
      file->write("transFluxRef", 1, &Nfreq, fluxTrans, false );
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

      minCrn = fluxBoxVol.get_min_corner();
      maxCrn = fluxBoxVol.get_max_corner();
      double crnFluxBox[nCrn] = {minCrn.x(),minCrn.y(),minCrn.z(),maxCrn.x(),maxCrn.y(),maxCrn.z()};
      file->write( "fluxBoxCrn", 1, &nCrn, crnFluxBox, false );
      file->prevent_deadlock();
    }
    else
    {
      file = f.open_h5file("data/sphereCrossSection", meep::h5file::READWRITE );
      double *fluxRef = refPlane.flux();
      int Nfreq = refPlane.Nfreq;
      file->write("scatPlaneFlux", 1, &Nfreq, fluxRef, false );
      file->prevent_deadlock();

      double *fluxTrans = transPlane.flux();
      file->write("transFluxScat", 1, &Nfreq, fluxTrans, false );
      file->prevent_deadlock();

      double *fluxOut = fluxBox.flux();
      file->write("scatFluxBox", 1, &Nfreq, fluxOut, false );
      file->prevent_deadlock();
      f.output_hdf5( meep::Dielectric, v.surroundings(), file, false, true );
      file->prevent_deadlock();
    }
  #else
    file = f.open_h5file("data/sphereCrossSection");
    f.output_hdf5( meep::Dielectric, v.surroundings(), file, false, true );
    file->prevent_deadlock();
    f.output_hdf5( fieldComp, v.surroundings(), file, false, true );
    file->prevent_deadlock();
  #endif
}
int main(int argc, char **argv) {
  meep::initialize mpi(argc, argv); // do this even for non-MPI Meep
  try
  {
    #ifdef GAUSSIAN_SOURCE
      run( epsRef );
      overwriteFile = false;
    #endif
    run( eps );
  }
  catch(...)
  {
    clog << "An exception occured...\n!";
  }
  return 0;
}
