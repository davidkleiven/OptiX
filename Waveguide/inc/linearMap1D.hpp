#ifndef LINEAR_MAP_1D_H
#define LINEAR_MAP_1D_H

/** Struct for handling points used in a linear map */
struct CorrespondingPoints
{
  double xFrom;
  double xTo;
};

/** Class that handles a linear map between two coordinate systems */
class LinearMap1D
{
public:
  LinearMap1D(){};

  /** Initialize the jacobian based on two coordinate sets that is known th match */
  void initialize( const CorrespondingPoints &p1, const CorrespondingPoints &p2 );

  /** Get the transformed coordinate */
  double get( double xFrom ) const;
private:
  double jacobian;
  double shift;
};
#endif
