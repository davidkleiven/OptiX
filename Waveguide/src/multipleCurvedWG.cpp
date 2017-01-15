#include "multipleCurvedWG.hpp"
#include <json/reader.h>
#include <stdexcept>
#include <fstream>
#include <cmath>
#include <cassert>
#include "curvedWGConfMap.hpp"
using namespace std;

const double PI = acos(-1.0);
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
  if ( solver != NULL ) delete solver;
  if ( src != NULL ) delete src;
}

void MultipleCurvedWG::loadWaveguides( const string &jsonfname )
{
  geometryfile = jsonfname;
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
  loadWaveguides( root );
}
void MultipleCurvedWG::loadWaveguides( const Json::Value &root )
{

  imagefile = root["figname"].asString();

  waveguides = new vector<CurvedWGConfMap*>;
  for ( unsigned int i=0;i<root["waveguides"].size();i++ )
  {
    double radius = root["waveguides"][i]["radius"].asDouble()*1E6;
    angles.push_back( root["waveguides"][i]["angle"].asDouble() );
    string curvature = root["waveguides"][i]["curvature"].asString();
    if ( i == 0 )
    {
      waveguides->push_back( new CurvedWGConfMap() );
    }
    else
    {
      CurvedWGConfMap *newwg = new CurvedWGConfMap();
      typedef CurvedWGConfMap::Curvature_t curv_t;
      if ( curvature == "convex" )
      {
        newwg->setCurvature( curv_t::CONVEX );
      }
      else
      {
        newwg->setCurvature( curv_t::CONCAVE );
      }
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

  solver = new CrankNicholson();
  PlaneWave* pw = new PlaneWave();
  double xmin = -2.0*params.at("width");
  double xmax = 2.0*params.at("width");
  double zmin = 0.0;
  double length = 0.0;
  const double PI = acos(-1.0);
  unsigned int totalNz = 0;
  cladding.setRefractiveIndex( params.at("delta"), params.at("beta") );

  for ( unsigned int i=0;i<waveguides->size();i++ )
  {
    length = (angles[i]*PI/180.0)*(*waveguides)[i]->getRadiusOfCurvature();
    (*waveguides)[i]->setWidth( params.at("width") );
    (*waveguides)[i]->setLongitudinalDiscretization( zmin, zmin+length, params.at("stepZ"), params.at("downSamplingZ") );
    (*waveguides)[i]->setTransverseDiscretization( xmin, xmax, params.at("stepX"), params.at("downSamplingX") );
    (*waveguides)[i]->setCladding( cladding );
    (*waveguides)[i]->setWaveLength( params.at("wavelength") );

    totalNz += (*waveguides)[i]->nodeNumberLongitudinal();
    zmin += length;
  }

  this->setTransverseDiscretization( xmin, xmax, params.at("stepX"), params.at("downSamplingX") );

  // NOTE: Now zmin is equal to the total length of the waveguide
  this->setLongitudinalDiscretization( 0.0, zmin, params.at("stepZ"), params.at("downSamplingZ") );
  this->setWaveLength( params.at("wavelength") );

  solver->setEquation( eq );
  solver->setBoundaryCondition( Solver2D::BC_t::TRANSPARENT );
  (*waveguides)[0]->setSolver( *solver );
  unsigned int Nx = (*waveguides)[0]->nodeNumberTransverse();

  transmittivity = new arma::vec( totalNz, arma::fill::zeros );
  intensity = new arma::mat( Nx/params.at("downSamplingX"), totalNz/params.at("downSamplingZ"), arma::fill::zeros );
  pw->setWavelength( params.at("wavelength") );
  src = pw;

  // Add some post processing module to the last waveguide
  farfield.setPadLength( pow(2,17) );
  farfield.setAngleRange( -1.0, 1.0 );
  *waveguides->back() << farfield << exitfield << exPhase;
}

void MultipleCurvedWG::solve()
{
  arma::cx_vec endSolution;
  arma::cx_vec *endSolutonAddress = &endSolution;

  fsource.setData( &endSolution ); // Stores a pointer to the endSolution array
  unsigned int counter = 0;
  cdouble im(0.0,1.0);
  for ( auto wg=waveguides->begin();wg != waveguides->end(); ++wg )
  {
    clog << "Running waveguide " << counter++ << endl;
    solver->reset();

    (*wg)->setSolver( *solver );
    if ( wg == waveguides->begin() )
    {
      (*wg)->setBoundaryConditions(*src );
    }
    else
    {
      double xmin = (*wg)->transverseDiscretization().min;
      double xmax = (*wg)->transverseDiscretization().max;

      if ( (*wg)->getCurvature() != (*(wg-1))->getCurvature() )
      {
        flipWrtCenterOfWG( endSolution );
        double phaseShift = phaseDifference( **(wg-1), **wg );
        endSolution *= exp( im*wavenumber*phaseShift );
      }

      (*wg)->setBoundaryConditions( fsource );
    }
    (*wg)->solve();
    processSolution( **wg );

    endSolution = solver->getLastSolution();

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
  unsigned int NcSol = solver->getSolution().n_cols;
  unsigned int nmax = NcInt > NcSol+NzNextFillStartIntensity ? NcSol:NcInt-NzNextFillStartIntensity;

  // Get a copy of the solution
  arma::mat intensitySolution =  abs( solver->getSolution() );
  if ( wg.getCurvature() == CurvedWGConfMap::Curvature_t::CONVEX )
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

  double intensityAtZero = (*waveguides)[0]->getTransmittivity().getIntensityAtZero();
  for ( unsigned int i=0;i<nmax;i++ )
  {
    double normalization = wg.getTransmittivity().getIntensityAtZero()/intensityAtZero;
    (*transmittivity)(i+NzNextFillStartTrans) = wg.getTransmittivity().get()[i]*normalization;
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

void MultipleCurvedWG::save( ControlFile &ctl )
{
  ParaxialSimulation::save( ctl );
  assert( maingroup != NULL );

  // Add radii and angles attributes
  hsize_t size = waveguides->size();
  H5::DataSpace attribSpace( 1, &size );
  H5::Attribute att = maingroup->createAttribute( "radius", H5::PredType::NATIVE_DOUBLE, attribSpace );

  vector<double> radii;
  for ( unsigned int i=0;i<waveguides->size();i++ )
  {
    radii.push_back( (*waveguides)[i]->getRadiusOfCurvature()/1E6 );
  }

  att.write( H5::PredType::NATIVE_DOUBLE, &radii[0] );

  att = maingroup->createAttribute( "angles", H5::PredType::NATIVE_DOUBLE, attribSpace );
  att.write( H5::PredType::NATIVE_DOUBLE, &angles[0] );

  // Add filenames read from the json geometry files
  H5::DataSpace stringDs( H5S_SCALAR );
  H5::StrType strdatatype( H5::PredType::C_S1, 256 );
  att = maingroup->createAttribute( "image", strdatatype, stringDs );
  att.write( strdatatype, imagefile );

  att = maingroup->createAttribute( "geofile", strdatatype, stringDs );
  att.write( strdatatype, geometryfile );

  // Store data
  saveArmaMat( *intensity, "amplitude", commonAttributes );
  saveArmaVec( *transmittivity, "transmittivity", commonAttributes );

  // Get the far field
  arma::vec res;
  farfield.result( waveguides->back()->getSolver(), res );
  vector<H5Attr> additionalAttrib;
  farfield.addAttrib( additionalAttrib );
  saveArmaVec( res, farfield.getName().c_str(), additionalAttrib );

  exitfield.result( waveguides->back()->getSolver(), res );
  saveArmaVec( res, exitfield.getName().c_str() );

  exPhase.result( waveguides->back()->getSolver(), res );
  saveArmaVec( res, exPhase.getName().c_str() );
}

double MultipleCurvedWG::phaseDifference( const CurvedWGConfMap &source, const CurvedWGConfMap &target ) const
{
  return 0.0; // Phase difference has no effect
  double r1 = source.getRadiusOfCurvature();
  double r2 = target.getRadiusOfCurvature();
  double distSource = source.longitudinalDiscretization().max;

  double phase = 0.0;
  if ( source.getCurvature() != target.getCurvature() )
  {
    phase = PI;
  }
  phase += (1.0 - r2/r1)*distSource;
  return phase;
}
