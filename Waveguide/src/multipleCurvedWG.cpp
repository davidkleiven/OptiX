#include "multipleCurvedWG.hpp"
#include "curvedWGConfMapGuideToGuide.hpp"
#include <json/reader.h>
#include <stdexcept>
#include <fstream>
#include <cmath>
#include <cassert>
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
  for ( unsigned int i=0;i<root["waveguides"].size();i++ )
  {
    double radius = root["waveguides"][i]["radius"].asDouble()*1E6;
    angles.push_back( root["waveguides"][i]["angle"].asDouble() );
    if ( i == 0 )
    {
      waveguides->push_back( new CurvedWGConfMap() );
    }
    else
    {
      int sign = radius > 0.0 ? 1:-1;
      CurvedWGConfMap *newwg = new CurvedWGConfMap();
      newwg->setSign( sign );
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
  double length = 0.0;
  const double PI = acos(-1.0);
  unsigned int totalNz = 0;
  cladding.setRefractiveIndex( params.at("delta"), params.at("beta") );

  for ( unsigned int i=0;i<waveguides->size();i++ )
  {
    if ( i == 0 ) zmin = 0.0;
    else zmin = length;
    length = (angles[i]*PI/180.0)*(*waveguides)[i]->getRadiusOfCurvature();
    (*waveguides)[i]->setWidth( params.at("width") );
    (*waveguides)[i]->setLongitudinalDiscretization( zmin, zmin+length, params.at("stepZ"), params.at("downSamplingZ") );
    (*waveguides)[i]->setTransverseDiscretization( xmin, xmax, params.at("stepX"), params.at("downSamplingX") );
    (*waveguides)[i]->setCladding( cladding );
    (*waveguides)[i]->setWaveLength( params.at("wavelength") );

    totalNz += (*waveguides)[i]->nodeNumberLongitudinal();
  }

  this->setTransverseDiscretization( xmin, xmax, params.at("stepX"), params.at("downSamplingX") );
  this->setLongitudinalDiscretization( 0.0, zmin+length, params.at("stepZ"), params.at("downSamplingZ") );

  solver.setEquation( eq );
  solver.setBoundaryCondition( Solver2D::BC_t::TRANSPARENT );
  (*waveguides)[0]->setSolver( solver );
  unsigned int Nx = (*waveguides)[0]->nodeNumberTransverse();

  transmittivity = new arma::vec( totalNz, arma::fill::zeros );
  intensity = new arma::mat( Nx/params.at("downSamplingX"), totalNz/params.at("downSamplingZ"), arma::fill::zeros );
  pw.setWavelength( params.at("wavelength") );
}

void MultipleCurvedWG::solve()
{
  arma::cx_vec endSolution;
  arma::cx_vec *endSolutonAddress = &endSolution;

  fsource.setData( &endSolution ); // Stores a pointer to the endSolution array
  unsigned int counter = 0;
  int prevSign = 1;
  for ( auto wg=waveguides->begin();wg != waveguides->end(); ++wg )
  {
    clog << "Running waveguide " << counter++ << endl;
    solver.reset();

    if ( (*wg)->getSign() != prevSign )
    {
      flipWrtCenterOfWG( endSolution );
    }

    (*wg)->setSolver( solver );
    if ( wg == waveguides->begin() )
    {
      (*wg)->setBoundaryConditions( pw );
    }
    else
    {
      double xmin = (*wg)->transverseDiscretization().min;
      double xmax = (*wg)->transverseDiscretization().max;
      fsource.setLimits( xmin, xmax );
      (*wg)->setBoundaryConditions( fsource );
    }
    (*wg)->solve();
    processSolution( **wg );

    endSolution = solver.getLastSolution();
    prevSign = (*wg)->getSign();

    // Just to be 100 % sure that the address does not change when reallocation is needed
    assert( endSolutonAddress == &endSolution );
  }

  // Set data to NULL pointer since the endSolution vector goes out of scope
  fsource.setData( NULL );
}

void MultipleCurvedWG::processSolution( CurvedWGConfMap &wg )
{
  assert( intensity != NULL );
  assert( transmittivity != NULL );

  // Set the upper limit
  unsigned int NcInt = intensity->n_cols;
  unsigned int NcSol = solver.getSolution().n_cols;
  unsigned int nmax = NcInt > NcSol+NzNextFillStartIntensity ? NcSol:NcInt-NzNextFillStartIntensity;

  // Get a copy of the solution
  arma::mat intensitySolution =  abs( solver.getSolution() );
  if ( wg.getSign() == -1 )
  {
    flipWrtCenterOfWG( intensitySolution );
  }

  // Copy intensity
  for ( unsigned int i=0;i<nmax;i++ )
  {
    for ( unsigned int j=0;j<intensity->n_rows;j++ )
    {
      (*intensity)(j,i+NzNextFillStartIntensity) = intensitySolution(j,i);
    }
  }

  NzNextFillStartIntensity += nmax;

  unsigned int thisTransSize = transmittivity->size();
  unsigned int computedTransSize = wg.getTransmittivity().get().size();
  nmax = thisTransSize > computedTransSize+NzNextFillStartTrans ? computedTransSize:thisTransSize-NzNextFillStartTrans;

  for ( unsigned int i=0;i<nmax;i++ )
  {
    (*transmittivity)(i+NzNextFillStartTrans) = wg.getTransmittivity().get()[i];
  }
  NzNextFillStartTrans += nmax;
}

template<class elemType>
void MultipleCurvedWG::flipWrtCenterOfWG( elemType array[], unsigned int N ) const
{
  double center = -(*waveguides)[0]->getWidth()/2.0;

  double xmin = xDisc->min;
  double xmax = xDisc->max;
  unsigned int ic = (center-xmin)*N/(xmax-xmin);
  unsigned int nflips = ic > N/2 ? N-ic:ic;
  for ( unsigned int i=1;i<nflips;i++ )
  {
    elemType copy = array[ic+i];
    array[ic+i] = array[ic-i];
    array[ic-i] = copy;
  }
}

void MultipleCurvedWG::flipWrtCenterOfWG( arma::cx_vec &vec ) const
{
  flipWrtCenterOfWG( vec.memptr(), vec.n_elem );
}

void MultipleCurvedWG::flipWrtCenterOfWG( arma::mat &mat ) const
{
  for ( unsigned int i=0;i<mat.n_cols;i++ )
  {
    flipWrtCenterOfWG( mat.colptr(i), mat.n_rows );
  }
}
