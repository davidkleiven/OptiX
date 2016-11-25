#include "exitFieldSource.hpp"
#include <complex>
#include <stdexcept>
#include <H5Cpp.h>

using namespace std;
typedef complex<double> cdouble;

void ExitFieldSource::load( H5::H5File &file )
{
  H5::DataSet amp = file.openDataSet("exitIntensity");
  H5::DataSet phase = file.openDataSet("exitPhase");
  H5T_class_t classType = amp.getTypeClass();
  H5::DataSpace dspace = amp.getSpace();
  hsize_t size = dspace.getSimpleExtentNpoints();
  arma::vec ampl(size);
  H5::DataSpace dMemspace(1,&size);
  amp.read( ampl.memptr(), H5::PredType::NATIVE_DOUBLE, dMemspace, dspace );

  dspace = phase.getSpace();
  arma::vec ph(size);
  phase.read( ph.memptr(), H5::PredType::NATIVE_DOUBLE, dMemspace, dspace );

  cdouble im(0.0,1.0);
  // % is elementwise multiplication
  values = ampl%arma::cos(ph) + im*ampl%arma::sin(ph);
  hasLoadedData = true;
  H5::Attribute attr = amp.openAttribute("xmin");
  H5::DataType type = attr.getDataType();
  attr.read( type, &xDisc.min );
  attr = amp.openAttribute("xmax");
  attr.read(type, &xDisc.max );
  xDisc.step = (xDisc.max-xDisc.min)/static_cast<double>(size);
}

void ExitFieldSource::setDiscretization( double min, double max )
{
  if ( !hasLoadedData )
  {
    throw (runtime_error("Need to load data first!"));
  }
  xDisc.min = min;
  xDisc.max = max;
  xDisc.step = (max-min)/values.n_elem;
}

unsigned int ExitFieldSource::closest( double x ) const
{
  return (x-xDisc.min)/xDisc.step;
}

cdouble ExitFieldSource::operator()( double x ) const
{
  if (( x >= xDisc.max-xDisc.step ) || ( x < xDisc.min ))
  {
    return 0.0;
  }
  unsigned int indx = closest(x);
  double x0 = xDisc.min + indx*xDisc.step;
  double x1 = x0 + xDisc.step;
  return ( (x1-x)*values(indx+1) + (x-x0)*values(indx) )/xDisc.step;
}
