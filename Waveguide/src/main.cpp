#include "waveGuideRadiusCurvature.hpp"
#include "cladding.hpp"
#include "numerovSolver.hpp"
#include <iostream>
#include <string>
#include <stdexcept>

using namespace std;

int main( int argc, char** argv )
{
  double eDensityTa = 4066.5; // nm^-3 same number as Salditt et al.
  Cladding cladding;
  cladding.setElectronDensity( eDensityTa );

  cout << "Potential in Ta: " << cladding.getPotential() << "nm^{-2}\n";

  WaveGuideLargeCurvature wg;
  wg.setRadiusOfCurvature( 40E6 ); // 40 mm
  wg.setWidth( 100.0 ); // 100 nm
  wg.setCladding( cladding );
  wg.setWaveLength( 0.1 );

  Numerov solver;
  solver.setUpperInitialCondition( 100.0, 1E-3, 0.0);
  solver.setLowerInitialCondition( -100.0, 0.0, 1E-3);
  solver.setStepsize(0.01); // 1 nm

  wg.setSolver( solver );
  try
  {
    solver.setPropgationWavenumberLimits( 0.999, 1.0);
    wg.solve();
    //string fname("data/waveguide");
    //wg.save(fname);
  }
  catch ( exception &exc )
  {
    cout << exc.what() << endl;
  }
  catch ( ... )
  {
    cout << "Unknown exception occured...\n";
  }
  return 0;
}
