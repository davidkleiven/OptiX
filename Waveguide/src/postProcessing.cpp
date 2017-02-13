#include "postProcessing.hpp"
#include "solver.hpp"
#include <cassert>

typedef post::PostProcessingModule::ReturnType_t ret_t;
typedef Solver::Dimension_t dim_t;
ret_t post::FieldQuantity::getReturnType( const Solver &solver ) const
{
  switch ( solver.getDim() )
  {
    case ( dim_t::TWO_D ):
      return ret_t::matrix2D;
    case ( dim_t::THREE_D ):
      return ret_t::cube3D;
  }
}

ret_t post::ProjectionQuantity::getReturnType( const Solver &solver ) const
{
  switch ( solver.getDim() )
  {
    case ( dim_t::TWO_D ):
      return ret_t::vector1D;
    case ( dim_t::THREE_D ):
      return ret_t::matrix2D;
  }
}

void post::PostProcessingModule::resizeMatrix( const arma::mat &orig, arma::mat &resized ) const
{
  assert( exportCols != 0 );
  assert( exportRows != 0 );
  if (( exportRows >= orig.n_rows ) || ( exportCols >= orig.n_cols ))
  {
    resized = orig;
    return;
  }

  double deltaX = static_cast<double>(orig.n_cols)/(static_cast<double>(exportCols)+1.0);
  double deltaY = static_cast<double>(orig.n_rows)/(static_cast<double>(exportRows)+1.0);
  resized.set_size( exportRows, exportCols );
  for ( unsigned int i=0;i<exportCols;i++ )
  {
    for ( unsigned int j=0;j<exportRows;j++ )
    {
      resized(j,i) = arma::mean( arma::mean(orig.submat(j*deltaY,i*deltaX, (j+1)*deltaY, (i+1)*deltaX)) );
    }
  }
}

void post::PostProcessingModule::setExportDimensions( unsigned int nrows, unsigned int ncols )
{
  exportRows = nrows;
  exportCols = ncols;
  resizeMatrices = true;
}
