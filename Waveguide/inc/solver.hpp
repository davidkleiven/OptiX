#ifndef SOLVER_BASE_CLASS_H
#define SOLVER_BASE_CLASS_H
#include <armadillo>
#include <visa/gaussianKernel.hpp>
#include <visa/lowPassFilter.hpp>

class ParaxialSimulation;

/** Base class for all solvers */
class Solver
{
public:
  Solver( const char* name ): name(name){};

  /** Get the name of the solution */
  std::string getName() const { return name; };

  /** Set waveguide to operate on */
  virtual void setSimulator( ParaxialSimulation &newGuide ){ guide = &newGuide; };

  /** Get the simulator */
  const ParaxialSimulation& getSimulator() const { return *guide; };

  /** Get solution */
  virtual const arma::cx_cube& getSolution3D() const;
  virtual const arma::cx_mat& getSolution() const;

  /** Propagates the solution one step */
  virtual void step() = 0;
protected:
  std::string name;
  ParaxialSimulation *guide{NULL};
  visa::GaussianKernel kernel;
  visa::LowPassFilter filter;
};

#endif
