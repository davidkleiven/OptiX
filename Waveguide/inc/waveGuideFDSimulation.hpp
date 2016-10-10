#ifndef WAVE_GUIDE_FD_SIMULATION_H
#define WAVE_GUIDE_FD_SIMULATION_H
#include <cstddef>
class Solver2D;

struct Disctretization
{
  double min;
  double max;
  double step;
};

class WaveGuideFDSimulation
{
public:
  WaveGuideFDSimulation();
  ~WaveGuideFDSimulation();
  void setTransverseDiscretization( double xmin , double xmax, double step );
  void setLongitudinalDiscretization( double zmin , double zmax, double step );
  unsigned int nodeNumberTransverse() const;
  unsigned int nodeNumberLongitudinal() const;
  const Disctretization& transverseDiscretization() const{ return *xDisc; };
  const Disctretization& longitudinalDiscretization() const { return *zDisc; };
private:
  Solver2D *solver{NULL};
  Disctretization *xDisc; // Transverse
  Disctretization *zDisc; // Along optical axis
};
#endif
