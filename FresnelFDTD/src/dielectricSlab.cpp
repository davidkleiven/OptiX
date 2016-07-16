#include "dielectricSlab.h"
#include <stdexcept>
#include <cmath>

const double PI = acos(-1.0);

const double DielectricSlab::xsize=5.0;
const double DielectricSlab::ysize=18.0;
const double DielectricSlab::pml_thick=3.0;
const double DielectricSlab::yc_plane=ysize/2.0;
const double DielectricSlab::source_y = ysize-pml_thick - 1.0;
double DielectricSlab::epslow = 1.0;
double DielectricSlab::epshigh = 1.0;

double DielectricSlab::kx=1.0;

DielectricSlab::DielectricSlab( double res ):vol(meep::vol2d(xsize, ysize,res)),sourcevol(NULL),struc(NULL),field(NULL){};

DielectricSlab::~DielectricSlab()
{
  if ( sourcevol != NULL ) delete sourcevol;
  if ( struc != NULL )     delete struc;
  if ( field != NULL )     delete field;
}
  

std::complex<double> DielectricSlab::amplitude( meep::vec &pos )
{
  const std::complex<double> IMAG_UNIT(0,1.0);
  return exp(IMAG_UNIT*DielectricSlab::kx*pos.x());
}

double DielectricSlab::dielectric( const meep::vec &pos )
{
  if ( pos.y() < yc_plane )
  {
    return epshigh;
  }
  return epslow;
}

void DielectricSlab::addSourceVol()
{ 
  meep::vec srcCorner1(0.0, source_y);
  meep::vec srcCorner2(xsize, source_y);
  if ( sourcevol != NULL ) delete sourcevol;
  sourcevol = new meep::volume(srcCorner1, srcCorner2);
}

void DielectricSlab::addStructure()
{
  if ( struc != NULL ) delete struc;
  struc = new meep::structure(vol, dielectric, meep::pml( pml_thick, meep::Y ));
}

void DielectricSlab::addField()
{
  if ( struc == NULL )
  {
    throw ( std::invalid_argument( "You must add a structure before adding fields!" ) );
  }
  if ( field != NULL ) delete field;
  field = new meep::fields( struc );
  field->use_bloch( meep::X, kx/(2.0*PI) ); // MEEP leaves out the factor 2pi (k = 1/lambda)
}

void DielectricSlab::setEpshigh( double epshigh )
{
  this->epshigh = epshigh;
}

void DielectricSlab::setEpslow( double epslow )
{
  this->epslow = epslow;
}

void DielectricSlab::setKx( double newKx )
{
  this->kx = newKx;
  if ( field != NULL )
  {
    field->use_bloch( meep::X, kx/(2.0*PI) );
  }
}

  const unsigned int NSTEPS = 20;

