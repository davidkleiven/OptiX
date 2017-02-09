
template<class wgPath, class wgShape>
void Waveguide3D::setShape( const wgShape &newshape )
{
  shape = newshape;
  shapeIsSet = true;
}

template<class wgPath, class wgShape>
void Waveguide3D::setPath( const wgPath &newPath )
{
  center = newPath;
  centerIsSet = true;
}

template<class wgPath, class wgShape>
bool Waveguide3D::isReady() const
{
  return shapeIsSet && centerIsSet;
}

template<class wgPath, class wgShape>
void Waveguide3D::setCladdingMaterial( const char* name )
{
  RefractiveIndex refr;
  refr.load(name);
  double energy = getEnergy();
  delta = refr.getDelta( energy );
  beta = refr.getBeta( energy );
}
