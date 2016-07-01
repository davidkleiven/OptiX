#ifndef FLUX_INTEGRATOR_H
#define FLUX_INTEGRATOR_H
#include <libscuff.h> 

class FluxIntegrator
{
  public:
    FluxIntegrator();
    ~FluxIntegrator();
    double compute( const scuff::RWGGeometry &geo, const IncField &IF, double omega, const double kBloch[2] );
    void setZpos(double zpos){_zpos=zpos;};
    void setnEvalPoints( unsigned int nPoints );
  private:
    double _zpos;
    unsigned int _nEvalPointsInEachDirection;
    HMatrix *_evaluationCrd;
    HMatrix *_flux;
};
    
#endif
