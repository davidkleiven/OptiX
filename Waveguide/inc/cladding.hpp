#ifndef CLADDING_MATERIAL_H
#define CLADDING_MATERIAL_H
#include <complex>
#include <vector>

typedef std::complex<double> cdouble;

/** Class handling the cladding of the waveguide */
class Cladding
{
public:
  Cladding(){};

  /** Set the electron density. Depricated */
  void setElectronDensity( double eDensity );

  /** Return the potential in the corresponding Schrodinger equation */
  double getPotential() const { return potential; };

  // Alternatives
  /** Set the refractive index */
  void setRefractiveIndex( double delta, double beta );

  /** Get delta. Refractive index n = 1 - delta + i*beta */
  double getDelta() const { return delta; }

  /** Get beta: Refractive index n=1 1- delta + i*beta */
  double getBeta() const { return beta; };
private:
  void computePotential();
  double delta{0.0};
  double beta{0.0};
  double thomsonScatteringLength{2.8179403267E-6}; // in nanometer
  double electronDensity{0.0};
  double potential{0.0}; // = 4pi*r0*rho
};
#endif
