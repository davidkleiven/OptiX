#include "absorber.hpp"

template<ApplyDim_t dim>
void Absorber::apply( arma::cx_mat &mat, DimGetter<dim> &getter ) const
{
  // Apply on both sides
  for ( unsigned int i=0;i<thickness;i++ )
  {
    getter(i,mat) *= exp( -static_cast<double>(thickness-i-1)*inverseDampingLength );
    getter(getter.size(mat)-i-1,mat) *= exp( -static_cast<double>(thickness-i-1)*inverseDampingLength );
  }
}

template void Absorber::apply<ApplyDim_t::COL>( arma::cx_mat &mat, DimGetter<ApplyDim_t::COL> &getter ) const;
template void Absorber::apply<ApplyDim_t::ROW>( arma::cx_mat &mat, DimGetter<ApplyDim_t::ROW> &getter ) const;
