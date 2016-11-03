#include <iostream>
#include "cladding.hpp"
#include "curvedWaveGuide2D.hpp"
#include "crankNicholson.hpp"
#include "controlFile.hpp"
#include "straightWG2D.hpp"
#include "planeWave.hpp"
#include "paraxialEquation.hpp"
#include "cylindricalParaxialEquation.hpp"
#include "curvedWGCylCrd.hpp"
#include "visualizer1D.hpp"
#include <complex>
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>

using namespace std;

int main( int argc, char **argv )
{
  bool cylindrical = false;
  for ( unsigned int i=1;i<argc;i++ )
  {
    string arg(argv[i]);
    if ( arg.find("--cyl") != string::npos )
    {
      cylindrical = true;
    }
  }

  srand( time(0) );
  const double PI = acos(-1.0);
  // Parameters for running a sweep over radii of curvature
  double LzOverR = 0.01; // max(z)/R << 1 is a requirement
  double xMarginAboveAndBelow = 0.5E3; // In nanometers = 0.5 um
  unsigned int Nz = 5000;
  unsigned int Nx = 5000;

  Cladding cladding;
  double delta = 4.49E-5;
  double beta = 3.45E-6;
  cladding.setRefractiveIndex(delta, beta);
  double width = 100.0; // Width of the waveguide in nm
  double wavelength = 0.1569;

  double Rcurv = 40E6;
  double zmin = 0.0;
  double zmax = Rcurv*LzOverR;
  double xmax = 2.0*width;
  double xmin = -1.0*width;
  //zmax = 600E3;

  CurvedWaveGuideFD wg;
  PlaneWave pw;
  ParaxialEquation eq;

  Visualizer visColor;
  Visualizer1D exVis;
  Visualizer1D accShift;
  CrankNicholson solver;

  try
  {
    visColor.init( "Top view" );
    exVis.init("Exit field");
    accShift.init("Accumulated Pixel Shift");
    solver.setEquation( eq );
    pw.setWavelength( wavelength );

    wg.setRadiusOfCurvature( Rcurv );
    wg.setWidth( width );
    wg.setWaveLength( wavelength );
    wg.setCladding( cladding );

    exVis.setLimits(-3.0,3.0);
    exVis.clear();

    double stepZ = (zmax-zmin)/static_cast<double>(Nz);
    double stepX = (xmax-xmin)/static_cast<double>(Nx);
    wg.setLongitudinalDiscretization(zmin,zmax,stepZ);
    clog << "N: " << Nz << " dx: " << stepX << " nm" << " dz: " << stepZ << " (rad or nm)\n";
    wg.setTransverseDiscretization(xmin,xmax,stepX);
    wg.setSolver(solver);
    wg.setBoundaryConditions( pw );
    wg.useBorderTracker();
    clog << "Border tracker initialized\n";
    wg.solve();
    clog << "Sytem solved\n";

    if ( exVis.isOpen() )
    {
      arma::vec exitF;
      wg.getExitField( exitF );
      exVis.fillVertexArray( exitF );
      exVis.display();
    }

    if ( accShift.isOpen() )
    {
      if ( wg.getBorderTracker() != NULL )
      {
        vector<int> vec = wg.getBorderTracker()->getAccumulatedPixelShift();
        arma::vec armVec( vec.size() );
        for ( unsigned int i=0;i<vec.size();i++ )
        {
          armVec(i) = vec[i];
        }
        double maxval = arma::max(armVec);
        double minval = arma::min(armVec);
        accShift.setLimits(minval,maxval);
        accShift.fillVertexArray( armVec );
        accShift.display();
      }
    }
    if ( visColor.isOpen() )
    {
      visColor.fillVertexArray( arma::abs( wg.getSolver().getSolution() ) );
      cout << "Maxval: " << arma::max( arma::abs( wg.getSolver().getSolution() ) ) << endl;
      visColor.display();
    }

    clog << " done\n";

    while ( visColor.isOpen() || exVis.isOpen() || accShift.isOpen() )
    {
      sf::Event event;
      while ( visColor.pollEvent( event ) )
      {
        if ( event.type == sf::Event::Closed )
        {
          visColor.close();
        }
      }

      while ( accShift.pollEvent( event ) )
      {
        if ( event.type == sf::Event::Closed )
        {
          accShift.close();
        }
      }

      while ( exVis.pollEvent( event ) )
      {
        if ( event.type == sf::Event::Closed )
        {
          exVis.close();
        }
      }
      exVis.display();
      visColor.display();
      accShift.display();
    }
  }
  catch ( exception &exc )
  {
    cout << exc.what() << endl;
    exVis.close();
    visColor.close();
    accShift.close();
    return 1;
  }

  return 0;
}
