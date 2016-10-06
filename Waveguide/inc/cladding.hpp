#ifndef CLADDING_MATERIAL_H
#define CLADDING_MATERIAL_H

class Cladding
{
public:
  Cladding();
  void setElectronDensity( double eDensity );
  void setThomsonPenetrationDepth( double depth );
  double getPotential() const { return potential; };
private:
  void computePotential();
  double delta{0.0};
  double beta{0.0};
  double thomsonPenetrationDepth{0.0};
  double electronDensity{0.0};
  double potential{0.0}; // = 4pi*r0*rho
};
#endif
