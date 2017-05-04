#include "sphere.hpp"
#include "planeWave.hpp"
#include "paraxialEquation.hpp"
#include "genericScattering.hpp"
#include "crankNicholson.hpp"
#include <iostream>
#include <visa/visa.hpp>
#include <pei/dialogBox.hpp>
#include <armadillo>
#include <chrono>
#include <thread>
#include <map>
#include <sstream>
#define KEEP_PLOT_FOR_SEC 3
//#define FAKE_RESULT_WITH_PURE_PHASE_SHIFT

using namespace std;

enum class Geometry_t{SPHERE, COATED_SPHERE};
int main( int argc, char **argv )
{
  cout << "--------------------\n";
  cout << "Modes:\n";
  cout << "0: Sphere\n";
  cout << "1: Coated sphere\n";
  cout << "--------------------\n";

  map<string,double> params;
  params["xmin_r"] = -4.0;            // Minimum x-position in units of the radius of the sphere
  params["xmax_r"] = 4.0;             // Maximum x-position in units of the radius of the sphere
  params["ymin_r"] = -4.0;            // Minimum y-position in units of the radius of the sphere
  params["ymax_r"] = 4.0;             // Maximum y-position in units of the radius of the sphere
  params["zmin_r"] = -2.0;            // Minimum z-position in units of the radius of the sphere
  params["zmax_r"] = 2.0;             // Maximum z-position in units of the radius of the sphere
  params["wavelength"] = 0.1569;      // Wavelength in nano meters
  params["radius"] = 2000.0;          // Radius of the sphere in nano meters
  params["Nt"] = 1024;                // Number of discretization points in the x and y direction (should 2^N)
  params["Nz"] = 512;                 // Number of discretization points in the z-direction (propagation direction)
  params["downSampleT"] = 8;          // Downsamplings factor in the x and y direction (1 --> no downsampling)
  params["downSampleZ"] = 4;          // Downsampling factor in the z-direction (1 --> no downsampling )
  params["savePlots"] = 0;            // Save the resulting plots to png file ( 0 --> no, otherwise yes )
  params["q_max_inverse_nm"] = 0.5;   // Maximum scattering wavevector to store for the farfield ( 1/nm )
  params["padLength"] = 2048;         // Length of the signal to use when computing the far field (should be 2^N)
  params["disableAbsorption"] = 0.0;  // Run without absorption (0 --> with absorption, otherwise run without absorption)
  params["mode"] = 0;
  params["coatThicknessIn_nm"] = 100.0; // Thickness of the coating in nano meters
  params["farFieldSaveSize"] = 512;     // Dimenstion of the farfield matrices to be saved
  params["intensityMinVis"] = 0.0;      // Minimum intensity for visualization
  params["intensityMaxVis"] = 1.1;      // Maximum intensity for visualization
  params["phaseMinVis"] = -1.58;        // Minimum phase for visualization
  params["phaseMaxVis"] = 1.58;         // Maximum phase for visualization

  pei::DialogBox dialog( params );
  dialog.show();

  Geometry_t geom;
  if ( static_cast<int>(params.at("mode")+0.5) == 0 )
  {
    geom = Geometry_t::SPHERE;
  }
  else if ( static_cast<int>(params.at("mode")+0.5) == 1 )
  {
    geom = Geometry_t::COATED_SPHERE;
  }

  double anglemax = params.at("q_max_inverse_nm")*params.at("wavelength")/(2.0*3.14159);
  anglemax *= (180.0/3.14159);

  // Define lengths in nm
  double r = params.at("radius");
  double xmin = params.at("xmin_r")*r;
  double xmax = params.at("xmax_r")*r;
  double zmin = params.at("zmin_r")*r;
  double zmax = params.at("zmax_r")*r;
  bool disableAbsorption = ( static_cast<int>(params.at("disableAbsorption")+0.5) != 0 );

  double dx = (xmax-xmin)/params.at("Nt");
  double dz = (zmax-zmin)/params.at("Nz");
  Point3D center;
  center.x = ( xmax+xmin)/2.0;
  center.y = center.x;
  center.z = (zmax+zmin)/2.0;

  GenericScattering simulation("sphere");
  stringstream ss;
  ss << "Simulation of a SiO2 sphere of radius " << r << " nm. ";
  simulation.setBeamWaist( 400.0*r );
  simulation.setMaxScatteringAngle( anglemax );
  simulation.xmin = xmin;
  simulation.xmax = xmax;
  simulation.ymin = xmin;
  simulation.ymax = xmax;
  simulation.zmin = zmin;
  simulation.zmax = zmax;
  simulation.dx = dx;
  simulation.dy = dx;
  simulation.dz = dz;
  simulation.downSampleX = params.at("downSampleT");
  simulation.downSampleY = params.at("downSampleT");
  simulation.downSampleZ = params.at("downSampleZ");
  simulation.wavelength = params.at("wavelength");
  simulation.phaseMin = params.at("phaseMinVis");
  simulation.phaseMax = params.at("phaseMaxVis");
  simulation.intensityMin = params.at("intensityMinVis");
  simulation.intensityMax = params.at("intensityMaxVis");
  simulation.realTimeVisualization = false;
  simulation.FFTPadLength = params.at("padLength");

  Sphere sphere( center, r );
  CoatedSphere coated( center, r, params.at("coatThicknessIn_nm") );

  if ( disableAbsorption )
  {
    ss << "Absorption was disabled. ";
    sphere.noAbsorption();
    coated.noAbsorption();
  }

  switch( geom )
  {
    case Geometry_t::SPHERE:
      simulation.setMaterial( sphere );
      break;
    case Geometry_t::COATED_SPHERE:
      ss << "The sphere was coated with a Au layer of thickness " << params.at("coatThicknessIn_nm") << " nm. ";
      simulation.setMaterial( coated );
      break;
  }

  try
  {
    sphere.setMaterial( "SiO2", simulation.getEnergy() );
    coated.setMaterial( "SiO2", simulation.getEnergy() );
    coated.setCoatingMaterial( "Au", simulation.getEnergy() );

    ss << "Delta sphere: " << sphere.getDelta() << ". Beta sphere: " << sphere.getBeta();
    if ( geom == Geometry_t::COATED_SPHERE )
    {
      ss << "Delta coating: " << coated.getDeltaCoating() << ". Beta coating: " << coated.getBetaCoating();
    }

    if ( argc > 1 )
    {
      simulation.imgname = argv[1];
    }

    simulation.description = ss.str();
    simulation.solve();
    simulation.save( "data/sphere.h5" );
  }
  catch ( exception &exc )
  {
    cout << exc.what() << endl;
    return 1;
  }
  catch (...)
  {
    cout << "Unrecognized exception!\n";
    return 1;
  }

  return 0;
}
