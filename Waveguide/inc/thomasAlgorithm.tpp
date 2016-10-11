
template <class T>
void ThomasAlgorithm::solve( T diag[], const T subdiag[], T rhs[], unsigned int N) const
{
  // Store a copy of the supdiag
  T supdiag[N-1];
  for ( unsigned int i=0;i<N-1;i++ )
  {
    supdiag[i] = subdiag[i];
  }

  // Forward step
  supdiag[0] /= diag[0];
  rhs[0] /= diag[0];
  for ( unsigned int i=1;i<N;i++ )
  {
    if ( i < N-1 )
    {
      supdiag[i] /= (diag[i] - subdiag[i-1]*supdiag[i-1]);
    }
    rhs[i] = (rhs[i] - subdiag[i-1]*rhs[i-1])/(diag[i] - subdiag[i-1]*supdiag[i-1]);
  }

  // Solve

  diag[N-1] = rhs[N-1];
  for ( int i=N-2;i>=0;i-- )
  {
    diag[i] = rhs[i] - supdiag[i]*diag[i+1];
  }
}
