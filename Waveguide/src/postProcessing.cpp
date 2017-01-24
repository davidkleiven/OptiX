#include "postProcessing.hpp"
#include "solver.hpp"

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
