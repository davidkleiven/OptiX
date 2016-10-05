#include "FormFactorFullSphere.h"
#include "GISASSimulation.h"
#include "Histogram2D.h"
#include "IntensityDataIOFactory.h"
#include "Layer.h"
#include "Materials.h"
#include "MultiLayer.h"
#include "Particle.h"
#include "ParticleLayout.h"
#include "Simulation.h"
#include "Units.h"
#include <string>
#include <sstream>

using namespace std;
enum class Mode_t {ONLY_SPHERE, SPHERE_ON_PLANE, SPHERE_ON_MEMBRANE};

int main(int argc, char **argv)
{
  string HELP_MSG("Usage: ./sphereSubstrate [arguments]\n");
  HELP_MSG += "--help: Print this message\n";
  HELP_MSG += "--sphere: Run simulation with a single sphere\n";
  HELP_MSG += "--onPlane: Run simulation with a sphere on a plane\n";
  HELP_MSG += "--onMembrane: Run simulation with a sphere on a membrane\n";
  HELP_MSG += "--thick: Thickness of the membrane in units of the radius of the sphere\n";

  // Parse arguments
  Mode_t mode = Mode_t::ONLY_SPHERE;
  double thicknessInUnitsOfRSphere = 0.1;
  for ( unsigned int i=1;i<argc;i++ )
  {
    string arg( argv[i] );
    if ( arg.find("--sphere") != string::npos )
    {
      mode = Mode_t::ONLY_SPHERE;
    }
    else if ( arg.find("--onPlane") != string::npos )
    {
      mode = Mode_t::SPHERE_ON_PLANE;
    }
    else if ( arg.find("--onMembrane") != string::npos )
    {
      mode = Mode_t::SPHERE_ON_MEMBRANE;
    }
    else if ( arg.find("--thick=") != string::npos )
    {
      stringstream ss;
      ss << arg.substr(7);
      ss >> thicknessInUnitsOfRSphere;
    }
    else if ( arg.find("--help") != string::npos )
    {
      cout << HELP_MSG;
      return 0;
    }
    else
    {
      cout << "Unknown argument " << arg << endl;
      return 0;
    }
  }

  // Refractive index of substrate
  // n = 1-delta + i*beta
  double refractiveIndexSubstrateDelta = 1.5E-5;
  double refractiveIndexSubstrateBeta = 2E-7;

  // Defina sample
  HomogeneousMaterial airMaterial("Air", 0.0, 0.0);
  HomogeneousMaterial substrateMaterial("Si3N4", refractiveIndexSubstrateDelta, refractiveIndexSubstrateBeta);

  // Define layers
  Layer airLayer(airMaterial);
  Layer substrateLayer(substrateMaterial);

  ParticleLayout partLayout;
  double refractiveIndexParticleDelta = 1E-5;
  double refractiveIndexParticleBeta = 1E-6;
  HomogeneousMaterial particleMaterial("Polymer", refractiveIndexParticleDelta, refractiveIndexParticleBeta);

  // Define a sphere with radius 1 um
  auto rSphere = 1.0*Units::micrometer;
  Particle sphere(particleMaterial, FormFactorFullSphere(rSphere));

  partLayout.addParticle(sphere);
  airLayer.addLayout(partLayout);

  MultiLayer sample;

  switch ( mode )
  {
    case Mode_t::ONLY_SPHERE:
      clog << "Running with only sphere in air\n";
      sample.addLayer(airLayer);
      break;
    case Mode_t::SPHERE_ON_PLANE:
      clog << "Running with a sphere on a plane\n";
      sample.addLayer(airLayer);
      sample.addLayer(substrateLayer);
      break;
    case Mode_t::SPHERE_ON_MEMBRANE:
      clog << "Running with a sphere on a membrane of finite thickness\n";
      sample.addLayer(airLayer);
      substrateLayer.setThickness(thicknessInUnitsOfRSphere*rSphere);
      sample.addLayer(substrateLayer);
      sample.addLayer(airLayer);
      break;
  }

  auto incidentGrazingAngle = 0.2*Units::degree;
  auto alpha_min = 0.15*Units::degree;
  auto alpha_max = 0.25*Units::degree;
  auto phi_min = -0.05*Units::degree;
  auto phi_max = 0.05*Units::degree;

  switch ( mode )
  {
    case Mode_t::ONLY_SPHERE:
      alpha_min = phi_min;
      alpha_max = phi_max;
      incidentGrazingAngle = 0.0*Units::degree;
      break;
  }

  // Setup Simulation
  GISASSimulation simulation;
  simulation.setDetectorParameters(1024, phi_min, phi_max, 1024, alpha_min, alpha_max);
  simulation.setBeamParameters(1.0*Units::angstrom, incidentGrazingAngle, 0.0*Units::degree);
  simulation.setSample(sample);

  // Run simulation
  simulation.runSimulation();
  const Histogram2D* result = simulation.getIntensityData();
  IntensityDataIOFactory::writeIntensityData(*result, "intensity.int");

  if ( result != NULL ) delete result;
  return 0;
}
