#ifndef WAVEGUIDE_3D_H
#define WAVEGUIDE_3D_H
#include "paraxialSimulation.hpp"
#include "refractiveIndex.hpp"

template<class wgPath, class wgShape>
class Waveguide3D: public ParaxialSimulation
{
public:
  Waveguide3D(): ParaxialSimulation("Waveguide3D"){};

  /** Sets the shape of the waveguide */
  void setShape( const wgShape &newShape );

  /** Set line at which the center of the waveguide follows */
  void setCenter( const wgPath &newPath );

  /** Set cladding material */
  void setCladdingMaterial( const char* name );
private:
    wgPath center;
    wgShape shape;
    double delta{0.0};
    double beta{0.0};

    bool centerIsSet{false};
    bool shapeIsSet{false};

    bool isReady() const;
};

// Include the implementation
#include "waveguide3D.tpp"
#endif
