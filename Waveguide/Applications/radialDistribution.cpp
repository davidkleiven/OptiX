#include <PaxPro/genericScattering.hpp>
#include <armadillo>
#include <iostream>
#include <cassert>

using namespace std;

class RadialDistribution: public MaterialFunction
{
public:
  RadialDistribution(){};

  void loadCSV( const char *fname )
  {
    arma::mat data;
    data.load( fname, arma::csv_ascii );
    radius = data.col(0);
    deltaArray = data.col(1);
    betaArray = data.col(2);
    rmin = radius.min();
    rmax = radius.max();
  }

  void getXrayMatProp( double x, double y, double z, double &delta, double &beta ) const override
  {
    double r = sqrt(x*x+y*y+z*z);

    int i1 =  (r-rmin)*(radius.n_elem-1)/( rmax-rmin );
    int i2 = i1+1;
    if ( i2 < radius.n_elem )
    {
      //cout << radius(i1) << " " << r << " " << radius(i2) << endl;
      //assert( r >= radius(i1) );
      //assert( r <= radius(i2) );
      double weight =  ( r-radius(i1) )/( radius(i2)-radius(i1) );

      if ( r <= radius(i1) ) weight = 0.0;
      if ( r >= radius(i2) ) weight = 1.0;
      assert( weight >= 0.0 );
      assert( weight <= 1.0 );

      delta = weight*deltaArray(i2) + (1.0-weight)*deltaArray(i1);
      beta = weight*betaArray(i2) + (1.0-weight)*betaArray(i1);
      assert( beta >= 0.0 );
    }
    else if ( i1==radius.n_elem-1 )
    {
      delta = deltaArray(i1);
      beta = betaArray(i1);
    }
    else
    {
      delta = 0.0;
      beta = 0.0;
    }
  }

  arma::vec radius;
  arma::vec deltaArray;
  arma::vec betaArray;

private:
  double rmin{0.0};
  double rmax{0.0};
};


// MAIN FUNCTION
int main( int argc, char** argv )
{
  if ( argc != 3 )
  {
    cout << "Usage: ./radialDistribution.out Nx Nz\n";
    return 3;
  }

  RadialDistribution material;
  material.loadCSV("data/radialSphereDist.csv");
  GenericScattering sim("RadialElectronDistribution");
  sim.description = "Simulation of a sphere consisting of a layered sphere 1) PMMA 2) Ni 3) Au 4) SiO2";
  sim.setBeamWaist( 400.0*material.radius.max() );

  sim.xmin = -1.1*material.radius.max();
  sim.xmax = 1.1*material.radius.max();
  sim.ymin = -1.1*material.radius.max();
  sim.ymax = 1.1*material.radius.max();
  sim.zmin = -1.1*material.radius.max();
  sim.zmax = 1.1*material.radius.max();

  unsigned int Nx = stoi( argv[1] );
  unsigned int Nz = stoi( argv[2] );

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
  double maxScatAngle = 0.5;
  sim.setMaxScatteringAngle( maxScatAngle );
  sim.wavelength = 0.177; // wavelength in nm

  sim.FFTPadLength = 8092;
  sim.setMaterial( material );
  unsigned int UID = rand()%1000000;
  sim.propagator = GenericScattering::SolverType_t::PROJ;
  try
  {
    sim.solve();

    stringstream fname;
    fname << "data/layeredSphere" << UID << ".h5";
    // Save results
    sim.save( fname.str() );
    clog << "Filename: " << fname.str() << endl;
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
