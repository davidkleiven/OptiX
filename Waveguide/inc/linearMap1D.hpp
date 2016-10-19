#ifndef LINEAR_MAP_1D_H
#define LINEAR_MAP_1D_H

struct CorrespondingPoints
{
  double xFrom;
  double xTo;
};

class LinearMap1D
{
public:
  LinearMap1D(){};
  void initialize( const CorrespondingPoints &p1, const CorrespondingPoints &p2 );
  double get( double xFrom ) const;
private:
  double jacobian;
  double shift;
};
#endif
