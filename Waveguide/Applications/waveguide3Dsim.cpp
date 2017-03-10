#include "waveguide3D.hpp"
#include "controlFile.hpp"
#include "postProcessMod.hpp"
#include "fftSolver3D.hpp"
#include "gaussianBeam.hpp"
#include "waveguideShapes.hpp"
#include "waveguide3DPaths.hpp"
#include <iostream>
#include <visa/visa.hpp>
#include <pei/dialogBox.hpp>
#include <map>
#include <string>
#include <SFML/Graphics.hpp>

using namespace std;
enum class Wg_t{ STRAIGHT, CURVED};

int main( int argc, char** argv )
{
  string imgStoring("");
  if ( argc > 1 )
  {
    imgStoring = argv[1];
  }
  cout << "--------------------\n";
  cout << "Modes:\n";
  cout << "0: Straight\n";
  cout << "1: Curved\n";
  cout << "--------------------\n";

  map<string,double> params;
  params["radiusInmm"] = 40.0;
  params["width"] = 100.0;
  params["depth"] = 100.0;
  params["Nt"] = 256;
  params["Nz"] = 64;
  params["wavelength"] = 0.1569;
  params["downSampleT"] = 8;
  params["downSampleZ"] = 4;
  params["q_max_inverse_nm"] = 0.5;
  params["padLength"] = 1024;
  params["zmax"] = 400E3;
  params["mode"] = 0;
  params["gaussianWaistInWidths"] = 1.0;
  params["yminInDepths"] = -2.0;
  params["ymaxInDepths"] = 5.0;
  params["absorberWidthInWavelengths"] = 20.0;
  params["absorberDampingInWavelengths"] = 4.0;

  pei::DialogBox box( params );
  box.show();

  double ymin = params.at("depth")*params.at("yminInDepths");
  double ymax = params.at("depth")*params.at("ymaxInDepths");
  double mode = params.at("mode");

  double anglemax = params.at("q_max_inverse_nm")*params.at("wavelength")/(2.0*3.14159);
  anglemax *= (180.0/3.14159);

  // Parse the mode
  Wg_t wgtype = Wg_t::STRAIGHT;
  if ( static_cast<int>(mode+0.5) == 0 )
  {
    wgtype = Wg_t::STRAIGHT;
  }
  else if ( static_cast<int>(mode+0.5) == 1 )
  {
    wgtype = Wg_t::CURVED;
  }
  else
  {
    cout << "Unknown mode!\n";
    return 1;
  }

  // Calculate the transverse boundaries
  double xmin, xmax;
  double width = params.at("width");

  switch ( wgtype )
  {
    case Wg_t::STRAIGHT:
      xmin = -6.0*width;
      xmax = 7.0*width;
      break;
    case Wg_t::CURVED:
      xmin = -pow( params.at("zmax"), 2)/(2.0*params.at("radiusInmm")*1E6 ) - width;
      xmax = 2.0*width;
      break;
  }

  // Post processing modules
  //post::Phase phase;
  //post::Intensity amplitude;
  post::ExitField ef;
  post::ExitIntensity ei;
  post::ExitPhase ep;
  post::FarField ff;

  // Set export size
  unsigned int exportRow = 512;
  unsigned int exportCol = 512;
  ef.setExportDimensions( exportRow, exportCol );
  ei.setExportDimensions( exportRow, exportCol );
  ff.setExportDimensions( exportRow, exportCol );

  SquareWell squareShape;
  squareShape.setWidth( params.at("width") );
  squareShape.setDepth( params.at("depth") );

  StraightPath strPath;
  ParabolicPath parPath;
  parPath.setRadius( params.at("radiusInmm") );

  Waveguide3D<StraightPath, SquareWell> straightWG;
  Waveguide3D<ParabolicPath, SquareWell> curvedWG;

  try
  {
    ff.setAngleRange( -anglemax, anglemax );
    ff.setPadLength( params.at("padLength")+0.5 );

    straightWG.setShape( squareShape );
    straightWG.setCenter( strPath );
    curvedWG.setShape( squareShape );
    curvedWG.setCenter( parPath );

    double dx = (xmax-xmin)/params.at("Nt");
    double dy = (ymax-ymin)/params.at("Nt");
    double dz = params.at("zmax")/params.at("Nz");


    FFTSolver3D solver;
    solver.overlayGeometry();

    if ( imgStoring != "" )
    {
      solver.storeImages( imgStoring.c_str() );
    }

    GaussianBeam gbeam;
    gbeam.setWavelength( params.at("wavelength") );
    gbeam.setDim( ParaxialSource::Dim_t::THREE_D );
    gbeam.setWaist( params.at("gaussianWaistInWidths")*params.at("width") );
    gbeam.setCenter( 0.0, 0.5*params.at("depth")/2.0 );

    solver.visualizeRealSpace();
    solver.setIntensityMinMax( 0.0, 2.0 );

    ControlFile ctl( "data/waveguide3D" );
    switch ( wgtype )
    {
      case Wg_t::STRAIGHT:
      {
        clog << "Running straight waveguide...\n";
        straightWG << ef << ei << ep << ff;
        straightWG.setWaveLength( params.at("wavelength") );
        straightWG.setTransverseDiscretization( xmin, xmax, dx, params.at("downSampleT") );
        straightWG.setVerticalDiscretization( ymin, ymax, dy );
        straightWG.setLongitudinalDiscretization( 0.0, params.at("zmax"), dz, params.at("downSampleZ") );
        straightWG.setCladdingMaterial( "MatProp/indexRefrTa.txt" );
        straightWG.setSolver( solver );

        double absorbWidth = params.at("absorberWidthInWavelengths")*params.at("wavelength");
        double damping = params.at("absorberDampingInWavelengths")*params.at("wavelength");
        solver.absorbingBC( absorbWidth, damping );

        straightWG.setBoundaryConditions( gbeam );
        straightWG.solve();
        straightWG.save( ctl );
        break;
      }
      case Wg_t::CURVED:
        clog << "Running curved waveguide...\n";
        curvedWG << ef << ei << ep << ff;
        curvedWG.setTransverseDiscretization( xmin, xmax, dx, params.at("downSampleT") );
        curvedWG.setVerticalDiscretization( ymin, ymax, dy );
        curvedWG.setLongitudinalDiscretization( 0.0, params.at("zmax"), dz, params.at("downSampleZ") );
        curvedWG.setCladdingMaterial( "MatProp/indexRefrTa.txt" );
        curvedWG.setSolver( solver );

        double absorbWidth = params.at("absorberWidthInWavelengths")*params.at("wavelength");
        double damping = params.at("absorberDampingInWavelengths")*params.at("wavelength");
        solver.absorbingBC( absorbWidth, damping );

        curvedWG.setBoundaryConditions( gbeam );
        curvedWG.solve();
        curvedWG.save( ctl );
        break;
    }
    ctl.save();
  }
  catch ( exception &exc )
  {
    cout << exc.what() << endl;
    return 0;
  }
  catch (...)
  {
    cout << "An unrecognized exception occured!\n";
    return 1;
  }


  return 0;
}
