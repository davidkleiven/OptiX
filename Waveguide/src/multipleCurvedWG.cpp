#include "multipleCurvedWG.hpp"
#include "curvedWGConfMapGuideToGuide.hpp"
#include <json/reader.h>
#include <stdexcept>
#include <fstream>
#include <cmath>
using namespace std;

MultipleCurvedWG::~MultipleCurvedWG()
{
  if ( waveguides != NULL )
  {
    for ( unsigned int i=0;i<waveguides->size();i++ )
    {
      delete (*waveguides)[i];
    }
    delete waveguides;
  }
  if ( intensity != NULL ) delete intensity;
  if ( transmittivity != NULL ) delete transmittivity;
}
void MultipleCurvedWG::loadWaveguides( const string &jsonfname )
{
  Json::Value root;
  Json::Reader reader;
  ifstream infile;
  infile.open( jsonfname.c_str() );
  if ( !infile.good() )
  {
    string msg("Could not open file ");
    msg += jsonfname;
    msg += "!";
    throw ( runtime_error(msg) );
  }
  reader.parse( infile, root );
  infile.close();

  waveguides = new vector<CurvedWGConfMap*>;

  Json::Value entries = root["waveguides"];
  double prevRadius;
  for ( unsigned int i=0;i<entries.size();i++ )
  {
    double radius = entries[i]["radius"].asDouble()*1E6;
    angles.push_back( entries[i]["angle"].asDouble() );
    if ( i == 0 )
    {
      waveguides->push_back( new CurvedWGConfMap() );
    }
    else
    {
      int sign = radius/(*waveguides)[0]->getRadiusOfCurvature() > 0.0 ? 1:-1;

      CurvedWGConfMapGuideToGuide *newwg = new CurvedWGConfMapGuideToGuide();
      newwg->setSign( sign );
      newwg->setTargetRadius( (*waveguides)[0]->getRadiusOfCurvature() );
      waveguides->push_back( newwg );
    }
    waveguides->back()->setRadiusOfCurvature( abs(radius) );
  }
}

void MultipleCurvedWG::init( const map<string,double> &params )
{
  if ( waveguides == NULL )
  {
    throw ( runtime_error("Waveguides needs to be loaded before the init function is called!") );
  }

  double xmin = -2.0*params.at("width");
  double xmax = 2.0*params.at("width");
  double zmin = 0.0;
  double zmax = 0.0;
  const double PI = acos(-1.0);
  unsigned int totalNz = 0;
  cladding.setRefractiveIndex( params.at("delta"), params.at("beta") );

  for ( unsigned int i=0;i<waveguides->size();i++ )
  {
    if ( i == 0 ) zmin = 0.0;
    else zmin = zmax;
    zmax = (angles[i]*PI/180.0)*(*waveguides)[i]->getRadiusOfCurvature();
    (*waveguides)[i]->setWidth( params.at("width") );
    (*waveguides)[i]->setLongitudinalDiscretization( zmin, zmax, params.at("stepZ"), params.at("downSamplingZ") );
    (*waveguides)[i]->setTransverseDiscretization( xmin, xmax, params.at("stepX"), params.at("downSamplingX") );
    (*waveguides)[i]->setCladding( cladding );

    totalNz += (*waveguides)[i]->nodeNumberLongitudinal();
  }
  solver.setEquation( eq );
  solver.setBoundaryCondition( Solver2D::BC_t::TRANSPARENT );
  (*waveguides)[0]->setSolver( solver );
  unsigned int Nx = (*waveguides)[0]->nodeNumberTransverse();

  transmittivity = new arma::vec( totalNz );
  intensity = new arma::mat( Nx/params.at("downSamplingX"), totalNz/params.at("downSamplingZ") );
  pw.setWavelength( params.at("wavelength") );
}

void MultipleCurvedWG::solve()
{
  arma::cx_vec endSolution;
  for ( auto wg=waveguides->begin();wg != waveguides->end(); ++wg )
  {
    (*wg)->setSolver( solver );
    if ( wg == waveguides->begin() )
    {
      (*wg)->setBoundaryConditions( pw );
    }
    else
    {
      solver.setLeftBC( endSolution.memptr() );
    }
    (*wg)->solve();
    processSolution( **wg );
    endSolution = solver.getSolution();
  }
}

void MultipleCurvedWG::processSolution( CurvedWGConfMap &wg )
{
  // Set the upper limit
  unsigned int NcInt = intensity->n_cols;
  unsigned int NcSol = solver.getSolution().n_cols;
  unsigned int nmax = NcInt > NcSol+NzNextFillStartIntensity ? NcSol:NcInt-NzNextFillStartIntensity;
  // Copy intensity
  for ( unsigned int i=0;i<nmax;i++ )
  {
    for ( unsigned int j=0;j<intensity->n_rows;j++ )
    {
      (*intensity)(j,i+NzNextFillStartIntensity) = abs( solver.getSolution()(j,i) );
    }
  }
  NzNextFillStartIntensity += nmax;

  unsigned int thisTransSize = transmittivity->size();
  unsigned int computedTransSize = wg.getTransmittivity().get().size();
  nmax = thisTransSize > computedTransSize+NzNextFillStartTrans ? computedTransSize:thisTransSize-NzNextFillStartTrans;

  for ( unsigned int i=0;i<nmax;i++ )
  {
    transmittivity[i+NzNextFillStartTrans] = wg.getTransmittivity().get()[i];
  }
  NzNextFillStartTrans += nmax;
}
