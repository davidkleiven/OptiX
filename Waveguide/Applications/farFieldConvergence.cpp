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
  bool useBorderTracker = true; // Does not seem to work properly yet
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
  const unsigned int NzDisc = 15;
  unsigned int Nz[NzDisc] = {250, 500, 1000, 2000, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000}; // Number of discretization points in x and z direction

  Cladding cladding;
  double delta = 4.49E-5;
  double beta = 3.45E-6;
  cladding.setRefractiveIndex(delta, beta);
  double width = 100.0; // Width of the waveguide in nm
  double wavelength = 0.1569;

  double Rcurv = 40E6;
  double zmin = 0.0;
  double zmax = Rcurv*LzOverR;
  double xmax = width+xMarginAboveAndBelow;
  double xmin = -0.5*zmax*LzOverR-xMarginAboveAndBelow;

  CurvedWaveGuideFD *wg = NULL;
  PlaneWave pw;
  ParaxialEquation *eq = NULL;

  if ( cylindrical )
  {
    wg = new CurvedWGCylCrd();
    CylindricalParaxialEquation *ceq = new CylindricalParaxialEquation();
    ceq->setRadiusOfCurvature( Rcurv );
    eq = ceq;
    zmax /= Rcurv;
    xmin = -width;
    xmax = 2.0*width;
  }
  else
  {
    wg = new CurvedWaveGuideFD();
    eq = new ParaxialEquation();
    if ( useBorderTracker )
    {
      xmin = -2.0*width;
      xmax = 3.0*width;
    }
  }

  Visualizer1D vis;
  Visualizer visColor;
  Visualizer1D exVis;
  CrankNicholson solver;

  try
  {
    vis.init( "Far field pattern" );
    visColor.init( "Top view" );
    exVis.init("Exit field");
    solver.setEquation( *eq );
    pw.setWavelength( wavelength );

    wg->setRadiusOfCurvature( Rcurv );
    wg->setWidth( width );
    wg->setWaveLength( wavelength );
    wg->setCladding( cladding );

    vis.setLimits(-6.0, 4.0);
    exVis.setLimits(-3.0,3.0);
    vis.clear();

    bool abortRun = false;
    for ( unsigned int iz=0;iz<NzDisc;iz++ )
    {
      double stepZ = (zmax-zmin)/static_cast<double>(Nz[iz]);
      double stepX = (xmax-xmin)/static_cast<double>(Nz[iz]);
      wg->setLongitudinalDiscretization(zmin,zmax,stepZ);
      clog << "N: " << Nz[iz] << " dx: " << stepX << " nm" << " dz: " << stepZ << " (rad or nm)\n";
      wg->setTransverseDiscretization(xmin,xmax,stepX);
      wg->setSolver(solver);
      wg->setBoundaryConditions( pw );
      wg->solve();
      wg->computeFarField( 65536 );
      if ( useBorderTracker )
      {
        wg->useBorderTracker();
      }

      arma::vec logFT = arma::log( wg->getFarField() );
      if ( iz == 0 )
      {
        double max = arma::max(logFT);
      //  vis.setLimits(max-8.0, max);
      }
      //arma::vec logFT = wg->getFarField();
      //cout << max( logFT ) << endl;

      double qMax = PI/stepX;
      double angleMax = qMax*180.0/(wg->getWavenumber()*PI);
      unsigned int indxStart = 0;
      unsigned int indxEnd = logFT.n_elem-1;
      if ( angleMax > 2.0 )
      {
        indxStart = (-2.0+angleMax)*static_cast<double>(logFT.n_elem-1)/(2.0*angleMax);
        indxEnd = (2.0+angleMax)*static_cast<double>(logFT.n_elem-1)/(2.0*angleMax);
        //indxStart = static_cast<double>(logFT.n_elem);
        //indxEnd = 0.6*static_cast<double>(logFT.n_elem);
        clog << "xmin: " << -2.0 << " " << "xmax: " << 2.0 << endl;
      }

      if ( vis.isOpen() )
      {
        vis.fillVertexArray( logFT.subvec( indxStart, indxEnd ) );
        vis.display();
      }

      if ( exVis.isOpen() )
      {
        arma::vec exitF;
        wg->getExitField( exitF );
        exVis.fillVertexArray( exitF );
        exVis.display();
      }

      if ( visColor.isOpen() )
      {
        visColor.fillVertexArray( arma::abs( wg->getSolver().getSolution() ) );
        visColor.display();
        sf::Event event;
        bool nextRun = false;
        clog << "Press space to continue\n";
        while (!nextRun && !abortRun)
        {
          while ( visColor.pollEvent( event ) )
          {
            if ( event.type == sf::Event::KeyPressed )
            {
              if ( event.key.code == sf::Keyboard::Space )
              {
                nextRun = true;
                clog << "Space pressed\n";
              }
              else if ( event.key.code == sf::Keyboard::Escape )
              {
                abortRun = true;
                vis.close();
                exVis.close();
                visColor.close();
              }
            }
          }
        }
      }
      clog << " done\n";
      if ( abortRun )
      {
        break;
      }
    }

    delete wg;
    delete eq;

    while ( vis.isOpen() || visColor.isOpen() || exVis.isOpen() )
    {
      sf::Event event;
      while ( vis.pollEvent( event ) )
      {
        if ( event.type == sf::Event::Closed )
        {
          vis.close();
        }
      }
      while ( visColor.pollEvent( event ) )
      {
        if ( event.type == sf::Event::Closed )
        {
          visColor.close();
        }
      }
      while ( exVis.pollEvent( event ) )
      {
        if ( event.type == sf::Event::Closed )
        {
          exVis.close();
        }
      }
      vis.display();
      exVis.display();
      visColor.display();
    }
  }
  catch ( exception &exc )
  {
    cout << exc.what() << endl;
    if ( vis.isOpen() )
    {
      vis.close();
    }
    delete wg;
    delete eq;
    return 1;
  }

  return 0;
}
