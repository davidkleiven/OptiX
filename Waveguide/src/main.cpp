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
#include "refractiveIndex.hpp"

using namespace std;

int main( int argc, char** argv )
{
  string coreMat("");
  string cladMat("SiO2");
  double energy = 10.0;
  for ( int i=1;i<argc;i++ )
  {
    string arg(argv[i]);
    if ( arg.find("--core=") != string::npos )
    {
      coreMat = arg.substr(7);
    }
    else if ( arg.find("--cladding=") != string::npos )
    {
      cladMat = arg.substr(11);
    }
    else if ( arg.find("--help") != string::npos )
    {
      cout << "Usage: ./waveguide --core=<corematerial> --cladding=<cladmat> --help\n";
      cout << "core: Material in the core (defaul:vacuum)\n";
      cout << "cladding: Material in the cladding (default=SiO2)\n";
      cout << "Energy: energy in keV (default=10keV)\n";
      return 1;
    }
    else if ( arg.find("--energy=") != string::npos )
    {
      stringstream ss;
      ss << arg.substr(9);
      ss >> energy;
    }
    else
    {
      cout << "Unknown command line argument " << arg << endl;
      return 1;
    }
  }
  RefractiveIndex refrIndex;
  RefractiveIndex refrCore;
  try
  {
    refrIndex.load(cladMat.c_str());
    if ( coreMat != "" )
    {
      refrCore.load(coreMat.c_str());
    }
  }
  catch ( exception &exc )
  {
    cout << exc.what() << endl;
    return 1;
  }
  gsl_error_handler_t* prevHandler = gsl_set_error_handler_off();
  double eDensityTa = 4066.5; // nm^-3 same number as Salditt et al.
  double beta = 3.45E-6; // Salditt is not use in the calculation, but useful when postprocessing
  double delta = 4.49E-5;
  double lambda = 0.1*12.398/energy;
  double width = 69.8;
  delta = refrIndex.getDelta( energy*1000.0 );
  beta = refrIndex.getBeta( energy*1000.0 );
  Cladding cladding;
  Cladding core;
  cladding.setElectronDensity( eDensityTa );
  cladding.setRefractiveIndex( delta, beta );

  cout << "Running width delta=" << cladding.getDelta() << endl;
  cout << "And beta = " << cladding.getBeta() << endl;

  //WaveGuideLargeCurvature wg;
  WaveGuideStraight wg;
  //wg.setRadiusOfCurvature( 40E6 ); // 40 mm
  wg.setWidth( width ); // 100 nm
  wg.setCladding( cladding );
  wg.setWaveLength( lambda );

  /** Set the material properties for the core */
  if ( coreMat != "" )
  {
    delta = refrCore.getDelta( energy*1000.0 );
    beta = refrCore.getBeta( energy*1000.0 );
    core.setRefractiveIndex( delta, beta );
    wg.setCore( core );
  }

  /*
  Numerov solver;
  solver.setUpperInitialCondition( 100, 1E-6, 0.0);
  solver.setLowerInitialCondition( -200.0, 0.0, 1E-6 );
  solver.setStepsize(0.0001); // 1 nm
  solver.setMaxIter(100000);
  */

  StandardFD solver;
  solver.setLimits( -2.0*width, width );
  solver.setStepsize(0.3);
  solver.setNumberOfModesToStore( 70 );

  wg.setSolver( solver );
  ControlFile ctl("data/eigenmodes");
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
