#ifndef VISUALIZER1D_H
#define VISUALIZER1D_H
#include "visualizer.hpp"

/** Class for visualizing 1D curves using the SFML library */
class Visualizer1D: public Visualizer
{
public:
  Visualizer1D(){};

  /** Set values to plot */
  void fillVertexArray( const arma::vec &vec );

  /** Set upper and lower limits on the y-axis */
  void setLimits( double minVal, double maxVal );
private:

  /** Get pixel corresponding to the y-value */
  unsigned int getY( double value ) const;
  double max{0.0};
  double min{0.0};
};
#endif
