#include "stokesParameters.hpp"
#include <stdexcept>

using namespace std;

unsigned int StokesParameters::componentIndx( meep::component c )
{
  if ( !firstComponentDetected && !secondComponentDetected )
  {
    comps[0] = c;
    firstComponentDetected = true;
    return 0;
  }
  else if ( firstComponentDetected && !secondComponentDetected )
  {
    if ( c == comps[0] ) return 0;
    else
    {
      comps[1] = c;
      secondComponentDetected = true;
      return 1;
    }
  }
  else
  {
    if ( c == comps[0] )
    {
      return 0;
    }
    else if ( c== comps[1] )
    {
      return 1;
    }
    else
    {
      throw ( runtime_error("Two components are detected and a third component is passed. Should not be the case for a flux plane!") );
    }
  }
}

void StokesParameters::computeIQ( meep::dft_chunk *Efield )
{
  I.resize(0);
  Q.resize(0);
  I.resize(Efield->Nomega, 0.0);
  Q.resize(Efield->Nomega, 0.0);
  unsigned int nElem = 0;

  // Compute I and Q
  for ( meep::dft_chunk* cur=Efield;cur!=NULL;cur=cur->next_in_dft )
  {
    int cIndx = componentIndx( cur->c );
    for ( unsigned int i=0;i<cur->N;i++ )
    {
      for ( unsigned int j=0;j<Efield->Nomega;j++ )
      {
        I[j] += abs(cur->dft[i*Efield->Nomega+j]);
        if ( cIndx == 0 )
        {
          Q[j] += abs(cur->dft[i*Efield->Nomega+j]);
        }
        else
        {
          Q[j] -= abs(cur->dft[i*Efield->Nomega+j]);
        }
      }
    }
    nElem += cur->N;
  }
  double Isum[Efield->Nomega];
  double Qsum[Efield->Nomega];

  // MPI stuff
  meep::sum_to_all( &I[0], Isum, Efield->Nomega );
  meep::sum_to_all( &Q[0], Qsum, Efield->Nomega );

  // Copy back and normalize
  for ( unsigned int i=0;i<Efield->Nomega;i++ )
  {
    I[i] = Isum[i]/(0.5*nElem);
    Q[i] = Qsum[i]/(0.5*nElem);
  }
}

void StokesParameters::computeUV( meep::dft_chunk* Efield )
{
  U.resize(0);
  V.resize(0);
  U.resize(Efield->Nomega, 0.0);
  V.resize(Efield->Nomega, 0.0);
  unsigned int nElem = 0;

  for ( meep::dft_chunk *cur=Efield; cur->next_in_dft != NULL; cur=cur->next_in_dft->next_in_dft )
  {
    int curCompIndx = componentIndx( cur->c );
    int nextCompIndx = componentIndx( cur->next_in_dft->c );
    if (( curCompIndx == nextCompIndx ) || (cur->N != cur->next_in_dft->N ))
    {
      throw( runtime_error("It is assumed that two adjacent chunks is the same volume, but different component!") );
    }
    for ( unsigned int i=0;i<cur->N;i++ )
    {
      for ( unsigned int j=0;j<cur->Nomega;j++ )
      {
        complex<double> val = cur->dft[i*cur->Nomega+j]*cur->next_in_dft->dft[i*cur->Nomega+j];
        U[j] += real(val);
        V[j] -= imag(val);
      }
    }
    nElem += cur->N;
  }

  double Usum[U.size()];
  double Vsum[V.size()];

  // MPI stuff
  meep::sum_to_all( &U[0], Usum, U.size() );
  meep::sum_to_all( &V[0], Vsum, V.size() );

  // Copy back land normalize
  for ( unsigned int i=0;i<U.size();i++ )
  {
    U[i] = Usum[i]/nElem;
    V[i] = Vsum[i]/nElem;
  }
}
