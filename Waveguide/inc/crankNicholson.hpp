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
  void solve() override final;
private:
  static unsigned int rowColToIndx( unsigned int Nz, unsigned int ix, unsigned int iz );
  static void indxToRowCol( unsigned int Nz, unsigned int indx, unsigned int &ix, unsigned int &iz );
  CSRBuild< std::complex<double> > matrix;

  unsigned int Nx;
  unsigned int Nz;
  double stepX, stepZ;
  double xmin;
  double zmin;
  double wavenumber;

  void solveCurrent( double iz );
};

#endif
