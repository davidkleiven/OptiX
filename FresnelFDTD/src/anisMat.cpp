#include "anisMat.h"
#include "dielectricSlab.h"

double StretchYMaterial::eps( const meep::vec &r )
{
  if ( slab->isInUpperHalfSpace( r ) )
  {
    if ( !isY(comp) )
    {
      return slab->getEpsUpper()/yscale;
    }
    else
    {
      return slab->getEpsUpper()*yscale;
    }
  }
  else
  {
    if ( !isY(comp) )
    {
      return slab->getEpsLower()/yscale;
    }
    else
    {
      return slab->getEpsLower()*yscale;
    }
  }
};

double StretchYMaterial::mu( const meep::vec &r )
{
  if ( !isY(comp) )
  {
    return 1.0/yscale;
  }
  else
  {
    return yscale;
  }
} 
