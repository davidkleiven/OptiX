#ifndef POST_PROCESSING_H
#define POST_PROCESSING_H
#include <string>
#include <armadillo>
#include "h5Attribute.hpp"

class Solver;
namespace post
{
class PostProcessingModule
{
public:
  enum class ReturnType_t { vector1D, matrix2D, cube3D };
  PostProcessingModule(const char* name ):name(name){};

  /** Returns the result */
  virtual void result( const Solver& solver, arma::cube& res ){};
  virtual void result( const Solver& solver, arma::mat& res ){};
  virtual void result( const Solver& solver, arma::vec& res ){};

  /** Add attribute to an external array */
  virtual void addAttrib( std::vector<H5Attr> &attrs ) const{};

  /** Returns the name */
  std::string getName() const { return name; };

  /** Sets the size of the exported matrices */
  void setExportDimensions( unsigned int rows, unsigned int cols );

  /** Which return type is used */
  virtual ReturnType_t getReturnType( const Solver& solver ) const = 0;
protected:
  std::string name;
  unsigned int exportRows{0};
  unsigned int exportCols{0};
  bool resizeMatrices{false};

  /** Resize the matrices to save space */
  void resizeMatrix( const arma::mat &mat, arma::mat &downsampled ) const;
};

/** A Field quantity is a quantity that is 2D for a 2D simulation and 3D for a 3D simulation
* Examples: Intensity, phase
* Counter examples: Average transmimittivity in a waveguide (1D in all cases), farField (1D for 2D sim and 2D for 3D sim)
*/
class FieldQuantity: public PostProcessingModule
{
public:
  FieldQuantity( const char* name ): PostProcessingModule(name){};
  virtual ReturnType_t getReturnType( const Solver& solver ) const override final;
};

/**
* A projection quantity is a quantity that has one lower dimension than the simulation
* Examples: far field
*/
class ProjectionQuantity: public PostProcessingModule
{
public:
  ProjectionQuantity( const char* name ): PostProcessingModule(name){};
  virtual ReturnType_t getReturnType( const Solver& solver ) const override final;
};
}; // namespace
#endif
