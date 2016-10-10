#ifndef CRANC_NICHOLSON_H
#define CRANC_NICHOLSON_H
#include "solver2D.hpp"
#include <complex>

template<class T>
struct CSRBuild
{
  T *elements{NULL};
  T *row{NULL};
  T *col{NULL};
};

class CrankNicholson: public Solver2D
{
public:
  CrankNicholson():Solver2D("CrankNicholson"){};
  ~CrankNicholson();
  void build() override final;
private:
  CSRBuild< std::complex<double> > matrix;
};

#endif
