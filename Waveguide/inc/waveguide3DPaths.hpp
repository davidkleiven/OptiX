#ifndef WAVEGUIDE3D_PATHS_H
#define WAVEGUIDE3D_PATHS_H

class StraightPath
{
public:
  /** Computes the x and y coordinate of the waveguide at position z */
  void get( double z, double &x, double &y ) const;

  /**
  This variable is included to be able to simulate a curved waveguide in cylindrical coordinates
  If R > 0.0, the quantity x/R is subtracted of the real part of the refractive index
  This is the first order correction to a straight waveguide, and is a good approximation
  when the radius of curvature is much larger than the width of the waveguide
  */
  double R{-10.0}; // Radius in millimeter
};

class ParabolicPath
{
public:
  /** Computes the x and y coordinate at position z */
  void get( double z, double &x, double &y ) const;

  /** Sets the radius of curvature */
  void setRadius( double rInmm ){ R = rInmm; };
private:
  double R{10.0}; // In mm
};
#endif
