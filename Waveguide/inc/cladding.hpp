#ifndef CLADDING_MATERIAL_H
#define CLADDING_MATERIAL_H

class Cladding
{
public:
  Cladding(){};
  void setElectronDensity( double eDensity );
  double getPotential() const { return potential; };
private:
  void computePotential();
  double delta{0.0};
  double beta{0.0};
  double thomsonScatteringLength{2.8179403267E-6}; // in nanometer
  double electronDensity{0.0};
  double potential{0.0}; // = 4pi*r0*rho
};
#endif
