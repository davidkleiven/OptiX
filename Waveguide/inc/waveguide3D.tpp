
#include <iostream>

template<class wgPath, class wgShape>
void Waveguide3D<wgPath,wgShape>::setShape( const wgShape &newshape )
{
  shape = newshape;
  shapeIsSet = true;
}

template<class wgPath, class wgShape>
void Waveguide3D<wgPath,wgShape>::setCenter( const wgPath &newPath )
{
  center = newPath;
  centerIsSet = true;
}

template<class wgPath, class wgShape>
bool Waveguide3D<wgPath,wgShape>::isReady() const
{
  return shapeIsSet && centerIsSet;
}

template<class wgPath, class wgShape>
void Waveguide3D<wgPath,wgShape>::setCladdingMaterial( const char* name )
{
  RefractiveIndex refr;
  std::string arg(name);

  // Quick fix to load file from a non standard location
  if ( arg.find("/") != std::string::npos )
  {
    // A file name was specified
    refr.loadUserDefinedFile( arg.c_str() );
  }
  else
  {
    refr.load(name);
  }
  double energy = getEnergy();
  delta = refr.getDelta( energy );
  beta = refr.getBeta( energy );
}

template<class wgPath, class wgShape>
void Waveguide3D<wgPath,wgShape>::getXrayMatProp( double x, double y, double z, double &matDelta, double &matBeta ) const
{
  double centerX;
  double centerY;
  center.get( z, centerX, centerY );

  if ( shape.isInside( x-centerX, y-centerY, z) )
  {
    matDelta = 0.0;
    matBeta = 0.0;
  }
  else
  {
    matDelta = delta;
    matBeta = beta;
  }
}
