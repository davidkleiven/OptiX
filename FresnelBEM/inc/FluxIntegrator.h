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
    void setYmin(double ymin){_ymin=ymin;};
    void setYmax(double ymax){_ymax=ymax;};
    void setXmin(double xmin){_xmin=xmin;};
    void setXmax(double xmax){_xmax=xmax;};
    void fillEvaluationXY();

    double incidentFlux( scuff::RWGGeometry &geo, IncField &IF, double omega, double kBloch[2] );
    double scatteredFlux( scuff::RWGGeometry &geo, HVector &vec, double omega, double kBloch[2] );
  private:
    double _zpos;
    unsigned int _nEvalPointsInEachDirection;
    double _ymin, _ymax, _xmin, _xmax;
    HMatrix *_evaluationCrd;
    HMatrix *_flux;
    double compute( scuff::RWGGeometry &geo, IncField *IF, HVector *vec, double omega, double kBloch[2] );
};
    
#endif
