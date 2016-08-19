#include "dielectricSlab.h"
#include <stdexcept>
#include <cmath>

const double PI = acos(-1.0);

const double DielectricSlab::xsize=5.0;
const double DielectricSlab::ysize=18.0;
const double DielectricSlab::pml_thick=3.0;
const double DielectricSlab::yc_plane=ysize/2.0;
const double DielectricSlab::source_y = ysize-pml_thick - 1.0;
double DielectricSlab::epsupper = 1.0;
double DielectricSlab::epslower = 1.0;
double DielectricSlab::kx=1.0;

DielectricSlab::DielectricSlab( double res ):vol(meep::vol2d(xsize, ysize,res)),sourcevol(NULL),struc(NULL),field(NULL), \
matX( meep::Ex, this ), matY( meep::Ey, this ), matZ( meep::Ez, this ){};

DielectricSlab::~DielectricSlab()
{
  if ( sourcevol != NULL ) delete sourcevol;
  if ( struc != NULL )     delete struc;
  if ( field != NULL )     delete field;
}
  

std::complex<double> DielectricSlab::amplitude( const meep::vec &pos )
{
  const std::complex<double> IMAG_UNIT(0,1.0);
  return exp(IMAG_UNIT*DielectricSlab::kx*pos.x());
}

double DielectricSlab::dielectric( const meep::vec &pos )
{
  if ( pos.y() < yc_plane )
  {
    return epslower;
  }
  return epsupper;
}

bool DielectricSlab::isInUpperHalfSpace( const meep::vec &pos )
{
  return pos.y() > yc_plane;
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
  bool anisAvg = false;
  double tol = 1.0; // --> set high since there is no averaging
  int maxeval = 1; // --> no averaging
  struc->set_chi1inv( meep::Ex, matX, false, tol, maxeval );
  struc->set_chi1inv( meep::Ey, matY, false, tol, maxeval );
  struc->set_chi1inv( meep::Ez, matZ, false, tol, maxeval );
  struc->set_chi1inv( meep::Hx, matX, false, tol, maxeval );
  struc->set_chi1inv( meep::Hy, matY, false, tol, maxeval );
  struc->set_chi1inv( meep::Hz, matZ, false, tol, maxeval );
  
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

void DielectricSlab::setEpsLower( double epslower )
{
  this->epslower = epslower;
}

void DielectricSlab::setEpsUpper( double epsupper )
{
  this->epsupper = epsupper;
}

void DielectricSlab::setKx( double newKx )
{
  this->kx = newKx;
  if ( field != NULL )
  {
    field->use_bloch( meep::X, kx/(2.0*PI) );
  }
}

void DielectricSlab::addSource( meep::src_time &source, meep::component fieldComp )
{
  if ( sourcevol == NULL )
  {
    throw ( std::invalid_argument("You have to add a souce volume before adding a source!" ) ); 
  }
  field->add_volume_source(fieldComp, source, *sourcevol, amplitude);
}

void DielectricSlab::output_hdf5( meep::component comp )
{
  if ( field == NULL ) return;
  
  field->output_hdf5( comp, vol.surroundings() );
}

void DielectricSlab::output_hdf5( meep::component comp, meep::h5file *file )
{
  if ( field == NULL ) return;
  field->output_hdf5( comp, vol.surroundings(), file );
} 
