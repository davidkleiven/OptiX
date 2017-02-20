#include "fieldMonitors.hpp"
#include <cmath>
#include <stdexcept>

using namespace std;

FieldMonitor::FieldMonitor( unsigned int N1, unsigned int N2 ): \
poynting( new arma::mat(N1, N2) ){};

FieldMonitor::~FieldMonitor()
{
  delete poynting;
}

void FieldMonitor::setDisplacementVectors( const meep::vec &v1, const meep::vec &v2 )
{
  if ( abs(v1&v2) > 1E-12 )
  {
    throw( runtime_error("Displacement vectors given to field monitor is not orthogonal" ) );
  }

  displacement1 = v1;
  displacement2 = v2;
}

void FieldMonitor::setIntensity( const meep::fields &field )
{
  for ( int i=0;i<poynting->n_cols;i++ )
  {
    for ( int j=0;j<poynting->n_rows;j++ )
    {
      meep::vec pos = origin + displacement1*(j-poynting->n_rows*0.5) + displacement2*(i-poynting->n_cols*0.5);
      (*poynting)(j,i) = 0.5*pow( abs(field.get_field( meep::Sx, pos )), 2 );
      (*poynting)(j,i) += 0.5*pow( abs(field.get_field( meep::Sy, pos )), 2 );
      (*poynting)(j,i) += 0.5*pow( abs(field.get_field( meep::Sz, pos )), 2 );
    }
  }
}
