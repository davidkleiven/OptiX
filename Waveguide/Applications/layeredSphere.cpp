#include <PaxPro/genericScattering.hpp>
#include <PaxPro/controlFile.hpp>
#include <iostream>
#include <sstream>
#include <ctime>

using namespace std;

class LayeredSphere: public MaterialFunction
{
public:
  void getXrayMatProp( double x, double y, double z, double &delta, double &beta ) const override
  {
    double r = sqrt( x*x+y*y+z*z );
    if ( r < rPolymer )
    {
      delta = deltaPoly;
      beta = betaPoly;
    }
    else if ( r < rNi )
    {
      delta = deltaNi;
      beta = betaNi;
    }
    else if ( r < rAu )
    {
      delta = deltaAu;
      beta = betaAu;
    }
    else if ( r < rSiO2 )
    {
      delta = deltaSiO2;
      beta = betaSiO2;
    }
    else
    {
      delta = 0.0;
      beta = 0.0;
    }
  };

  double getRmax() const { return rSiO2; };
private:
  double rPolymer{1541.0}; // nm
  double rNi{1597.0};      // nm
  double rAu{1630.0};      // nm
  double rSiO2{1691.0};    // nm
  double deltaPoly{5.46e-6};
  double betaPoly{1.59e-8};
  double deltaNi{3.39e-5};
  double betaNi{8.71e-7};
  double deltaAu{6.27e-5};
  double betaAu{8.04e-6};
  double deltaSiO2{1.48e-5};
  double betaSiO2{2.50e-7};
};

int main( int argc, char **argv )
{
  srand(time(0));
  if ( argc != 4 )
  {
    cout << "Usage: ./layeredSphere.out <fft or adim> <Nx> <Nz>\n";
    return 0;
  }
  LayeredSphere material;

  GenericScattering sim("LayeredSphere");
  sim.description = "Simulation of a sphere consisting of a layered sphere 1) PMMA 2) Ni 3) Au 4) SiO2";
  sim.setBeamWaist( 400.0*material.getRmax() );
  sim.xmin = -1.1*material.getRmax();
  sim.xmax = 1.1*material.getRmax();
  sim.ymin = -1.1*material.getRmax();
  sim.ymax = 1.1*material.getRmax();
  sim.zmin = -1.1*material.getRmax();
  sim.zmax = 1.1*material.getRmax();
  unsigned int Nx = stoi( argv[2] );
  unsigned int Nz = stoi( argv[3] );
  sim.dx = (sim.xmax-sim.xmin)/Nx;
  sim.dy = (sim.ymax-sim.ymin)/Nx;
  sim.dz = (sim.zmax-sim.zmin)/Nz;

  // Set size of exported files
  sim.exportNx = 2048;
  sim.exportNy = 2048;

  unsigned int finalSizeX = 256;
  unsigned int finalSizeZ = 256;
  sim.downSampleX = Nx/finalSizeX;
  sim.downSampleY = sim.downSampleX;
  sim.downSampleZ = Nz/finalSizeZ;
  sim.setMaxScatteringAngle( 0.5 );
  sim.wavelength = 0.177; // wavelength in nm

  string solver(argv[1]);
  sim.useFFTSolver = ( solver == "fft" );

  sim.FFTPadLength = 32768;
  sim.setMaterial( material );

  try
  {
    sim.solve();

    stringstream fname;
    fname << "data/layeredSphere" << rand()%1000000 << ".h5";
    // Save results
    sim.save( fname.str() );
  }
  catch( exception &exc )
  {
    cout << exc.what() << endl;
    return 1;
  }
  catch( ... )
  {
    cout << "Unrecognized exception...\n";
    return 2;
  }
  return 0;
}
