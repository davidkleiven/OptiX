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

int main(int argc, char **argv)
{
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
  Particle sphere(particleMaterial, FormFactorFullSphere(1.0*Units::micrometer));

  partLayout.addParticle(sphere);
  airLayer.addLayout(partLayout);

  MultiLayer sample;
  sample.addLayer(airLayer);
  sample.addLayer(substrateLayer);

  // Setup Simulation
  GISASSimulation simulation;
  simulation.setDetectorParameters(400, -1.0*Units::degree, 1.0*Units::degree, 400, 0.0*Units::degree, 2.0*Units::degree);
  simulation.setBeamParameters(1.0*Units::angstrom, 0.2*Units::degree, 0.0*Units::degree);
  simulation.setSample(sample);

  // Run simulation
  simulation.runSimulation();
  const Histogram2D* result = simulation.getIntensityData();
  IntensityDataIOFactory::writeIntensityData(*result, "intensity.int");

  if ( result != NULL ) delete result;
  return 0;
}
