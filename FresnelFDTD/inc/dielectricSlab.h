#ifndef DIELECTRIC_SLAB_H
#define DIELECTRIC_SLAB_H
#include "meep.hpp"

/**
* @class DielectricSlab  abstracts the underlying geometry parameters of the slab geometry.
*/
class DielectricSlab
{
public:
  explicit DielectricSlab(double resolution);
  ~DielectricSlab();
  inline meep::volume& getSourceVol()    {return *sourcevol;};
  inline meep::structure& getStructure() {return *struc;};
  inline meep::fields& getField()       {return *field;};
  void addSourceVol();
  void addStructure();
  void addField();
  void setKx(double kx);
  void setEpshigh(double epshigh);
  void setEpslow(double epslow);
  void addSource( meep::src_time &source, meep::component fieldComp );
  inline double getXsize() const {return xsize;};
  inline double getYsize() const {return ysize;};
  inline double getPMLThickness() const {return pml_thick;};
  inline double getSourceY() const {return source_y;};
  inline double getYcPlane() const {return yc_plane;};
  inline double getEpsLow() const {return epslow;};
  inline double getEpsHigh() const {return epshigh;};
  void output_hdf5( meep::component comp );
  void output_hdf5( meep::component comt, meep::h5file* file ); 
private:
  meep::grid_volume vol;
  meep::volume *sourcevol;
  meep::structure *struc;  
  meep::fields *field;
  static const double xsize;
  static const double ysize;
  static const double pml_thick;
  static const double yc_plane;
  static const double source_y;
  static double kx;
  static double epslow;
  static double epshigh;
  static std::complex<double> amplitude( const meep::vec &pos );
  static double dielectric( const meep::vec &pos ); 
};  
#endif
