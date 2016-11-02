#ifndef VISUALIZER1D_H
#define VISUALIZER1D_H
#include "visualizer.hpp"

class Visualizer1D: public Visualizer
{
public:
  Visualizer1D(){};
  void fillVertexArray( const arma::vec &vec );
  void setLimits( double minVal, double maxVal );
private:
  unsigned int getY( double value ) const;
  double max{0.0};
  double min{0.0};
};
#endif
