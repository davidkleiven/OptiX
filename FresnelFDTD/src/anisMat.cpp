#include "anisMat.h"
#include "dielectricSlab.h"
#include <iostream>

double StretchYMaterial::eps( const meep::vec &r )
{
  static bool printMsg = true;
  if ( printMsg )
  {
    std::cout << "Was in eps function...\n";
    printMsg = false;
  } 
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
  static bool printMsg = true;
  if ( printMsg )
  {
    std::cout << "Was in mu function...\n";
    printMsg = false;
  } 
  if ( !isY(comp) )
  {
    return 1.0/yscale;
  }
  else
  {
    return yscale;
  }
} 

void StretchYMaterial::eff_chi1inv_row( meep::component c, double chi1inv_row[3], const meep::volume &v, double tol, \
    int maxeval)
{
  meep::vec pos = (v.get_min_corner()+v.get_max_corner())*0.5;
  double epsilon = slab->dielectric(pos);
  double mu = 1.0;
  switch (c)
  {
    case meep::Ex:
      chi1inv_row[0] = 1.0/(epsilon/yscale);
      chi1inv_row[1] = 0.0;
      chi1inv_row[2] = 0.0;
      return;
    case meep::Ey:
      chi1inv_row[0] = 0.0;
      chi1inv_row[1] = 1.0/(epsilon*yscale);
      chi1inv_row[2] = 0.0;
      return;
    case meep::Ez:
      chi1inv_row[0] = 0.0;
      chi1inv_row[1] = 0.0;
      chi1inv_row[2] = 1.0/(epsilon/yscale);
      return;
    case meep::Hx:
      chi1inv_row[0] = 1.0/(mu/yscale);
      chi1inv_row[1] = 0.0;
      chi1inv_row[2] = 0.0;
      return;
    case meep::Hy:
      chi1inv_row[0] = 0.0;
      chi1inv_row[1] = 1.0/(mu*yscale);
      chi1inv_row[2] = 0.0;
      return;
    case meep::Hz:
      chi1inv_row[0] = 0.0;
      chi1inv_row[1] = 0.0;
      chi1inv_row[2] = 1.0/(mu/yscale);
    default:
      return;
  }
}
