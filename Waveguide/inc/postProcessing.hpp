#ifndef POST_PROCESSING_H
#define POST_PROCESSING_H
#include <string>
#include <armadillo>

class Solver2D;
namespace post
{
class PostProcessingModule
{
public:
  PostProcessingModule(const char* name):name(name){};

  /** Returns the result */
  virtual void result( const Solver2D& solver, arma::mat& res ){};
  virtual void result( const Solver2D& solver, arma::vec& res ){};

  /** Returns the name */
  std::string getName() const { return name; };
protected:
  std::string name;
};
};
#endif
