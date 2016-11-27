#ifndef CYLINDRICAL_PARAXIAL_EQUATION_H
#define CYLINDRICAL_PARAXIAL_EQUATION_H
#include "paraxialEquation.hpp"

/** Class for handling the coefficient in the cylindrical paraxial wave equation */
class CylindricalParaxialEquation: public ParaxialEquation
{
public:
  CylindricalParaxialEquation(){};

  /** Return the coefficient in front of first derivative */
  double F( double r, double theta ) const override final;

  /** Returns the coefficient in front of the the Laplacian */
  double G( double r, double theta ) const override final;

  /** Returns the coefficient inside first derivative in the Laplacian */
  double H( double r, double theta ) const override final;

  /** Returns extra refractive index that is added to delta */
  double J( double r, double theta ) const override final;

  /** Sets the radius of curvature in nano meter */
  void setRadiusOfCurvature( double Rcurv ){ R = Rcurv; };

  /** Returns the phase factor at a location theta. Depricated */
  virtual cdouble phaseFactor( double k, double theta ) const override;
private:
    double R{1.0};
};

#endif
