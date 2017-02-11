#ifndef WAVEGUIDE_SHAPES_H
#define WAVEGUIDE_SHAPES_H

/**
@brief: Describes a square well. The waveguide occupies the region:
y > 0.0 and x in (-width/2.0,width/2.0).
The region y > depth also belongs to the waveguide
*/
class SquareWell
{
public:
  SquareWell(){};

  /** Returns true if the coordinate belonds to the shape */
  bool isInside( double x, double y, double z ) const;

  /** Set the width of the waveguide in nm */
  void setWidth( double widthInm ){ width = widthInm; };

  /** Set the depth of the waveguide in nm */
  void setDepth( double depthInm ){ depth = depthInm; };
private:
  double width{1.0}; // In nm
  double depth{1.0};
};
#endif
