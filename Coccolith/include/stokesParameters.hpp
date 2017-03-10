#ifndef STOKES_PARAMETERS_H
#define STOKES_PARAMETERS_H
#include <meep.hpp>
#include <vector>

class StokesParameters
{
public:
  StokesParameters(){};

  /** Computes the Fourier Transformed stokes parameters from the fields */
  void compute( meep::dft_chunk *Efield );

  /** Get the Fourier transformed I */
  std::vector<double>& getI() { return I; };

  /** Get the Fourier transformed Q */
  std::vector<double>& getQ() { return Q; };

  /** Get the Fourier transformed U */
  std::vector<double>& getU() { return U; };

  /** Get the Fourier transformed V */
  std::vector<double>& getV() { return V; };
private:
  meep::component comps[2];

  bool firstComponentDetected{false};
  bool secondComponentDetected{false};

  // Stokes parameters see definition on wikipedia: https://en.wikipedia.org/wiki/Stokes_parameters
  std::vector<double> I;
  std::vector<double> Q;
  std::vector<double> U;
  std::vector<double> V;

  /** Returns the index of the present component */
  unsigned int componentIndx( meep::component c );

  /** Compute the two first Stokes parameters */
  void computeIQ( meep::dft_chunk *Efield );

  /** Computes the two last Stokes parameters */
  void computeUV( meep::dft_chunk *Efield );
};
#endif
