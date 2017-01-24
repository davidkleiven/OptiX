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
  enum class ReturnType_t { vector1D, matrix2D };
  PostProcessingModule(const char* name, ReturnType_t type ):name(name), returnType(type){};

  /** Returns the result */
  virtual void result( const Solver& solver, arma::cube& res ){};
  virtual void result( const Solver& solver, arma::mat& res ){};
  virtual void result( const Solver& solver, arma::vec& res ){};

  /** Add attribute to an external array */
  virtual void addAttrib( std::vector<H5Attr> &attrs ) const{};

  /** Returns the name */
  std::string getName() const { return name; };

  /** Which return type is used */
  ReturnType_t getReturnType() const { return returnType; };
protected:
  std::string name;
  ReturnType_t returnType;
};
};
#endif
