#ifndef WAVEGUIDE3D_PATHS_H
#define WAVEGUIDE3D_PATHS_H

class StraightPath
{
public:
  /** Computes the x and y coordinate of the waveguide at position z */
  void get( double z, double &x, double &y ) const;
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
