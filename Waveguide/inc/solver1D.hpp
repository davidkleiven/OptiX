#ifndef SOLVER_H
#define SOLVER_H
#include <string>

class Solver1D
{
  public:
    Solver1D( const char* name ): name(name){};
    virtual void solve() = 0;
  protected:
    std::string name;
};
#endif
