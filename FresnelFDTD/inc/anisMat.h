#ifndef ANIS_MAT_H
#define ANIS_MAT_H
#include "meep.hpp"

//#include "dielectricSlab.h"
class DielectricSlab;

class StretchYMaterial: public meep::material_function
{
  public:
    StretchYMaterial( const DielectricSlab *slab ): yscale(1.0), comp(comp), slab(slab){};
    inline double getYscale() const { return yscale; };
    void setYscale( double newyscale ){ yscale = newyscale; };
    
    // Override eps
    virtual double eps( const meep::vec &r );
    virtual double mu( const meep::vec &r );
    virtual bool has_mu() { return true; };
    virtual void eff_chi1inv_row( meep::component c, double chi1inv_row[3], const meep::volume &v, double tol, int maxeval);
  private:
    double yscale;
    meep::component comp;
    const DielectricSlab *slab;

    // Helper
    bool isX( meep::component comp) { return (comp==meep::Ex) || (comp==meep::Hx) || (comp==meep::Dx) || (comp==meep::Bx); };
    bool isY( meep::component comp) { return (comp==meep::Ey) || (comp==meep::Hy) || (comp==meep::Dy) || (comp==meep::By); };
    bool isZ( meep::component comp) { return (comp==meep::Ez) || (comp==meep::Hz) || (comp==meep::Dz) || (comp==meep::Bz); };
};

#endif
