#include "waveGuideRadiusCurvature.hpp"
#include "cladding.hpp"
#include "numerovSolver.hpp"
#include "waveGuideStraight.hpp"
#include <iostream>
#include <string>
#include <stdexcept>
#include <gsl/gsl_errno.h>
#include "stdFDsolver.hpp"
#include "controlFile.hpp"

using namespace std;

int main( int argc, char** argv )
{
  gsl_error_handler_t* prevHandler = gsl_set_error_handler_off();
  double eDensityTa = 4066.5; // nm^-3 same number as Salditt et al.
  Cladding cladding;
  cladding.setElectronDensity( eDensityTa );

  cout << "Potential in Ta: " << cladding.getPotential() << "nm^{-2}\n";

  WaveGuideLargeCurvature wg;
  //WaveGuideStraight wg;
  wg.setRadiusOfCurvature( 40E6 ); // 40 mm
  wg.setWidth( 100.0 ); // 100 nm
  wg.setCladding( cladding );
  wg.setWaveLength( 0.1569 );

  /*
  Numerov solver;
  solver.setUpperInitialCondition( 100, 1E-6, 0.0);
  solver.setLowerInitialCondition( -200.0, 0.0, 1E-6 );
  solver.setStepsize(0.0001); // 1 nm
  solver.setMaxIter(100000);
  */

  StandardFD solver;
  solver.setLimits( -150.0, 50.0 );
  solver.setStepsize(1.0);
  solver.setNumberOfModesToStore( 10 );

  wg.setSolver( solver );
  ControlFile ctl;
  try
  {
    //solver.setPropgationWavenumberLimits( 0.0, 0.5*cladding.getPotential() );
    wg.solve();
    wg.save(ctl);
    ctl.save();
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
