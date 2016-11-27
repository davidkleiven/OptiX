#ifndef THOMAS_ALGORITHM_H
#define THOMAS_ALGORITHM_H

/** Class for solving tridiagonal systems of equations */
class ThomasAlgorithm
{
public:
  ThomasAlgorithm(){};

  // On return the solution is stored in diag. NOTE: rhs is modified
  /** Solve the system of equations */
  template <class T>
  void solve( T diag[], const T subdiag[], T rhs[], unsigned int N) const;
};

#include "thomasAlgorithm.tpp"
#endif
