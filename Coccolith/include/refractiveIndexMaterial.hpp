#ifndef REFRACTIVE_INDEX_MATERIAL_H
#define REFRACTIVE_INDEX_MATERIAL_H
#include <vector>
#include "config.h"
#include <json/reader.h>

struct Lorentzian
{
  double preFactor{0.0};
  double resonanceWavelength{1.0};
};

/** This is a class that have dispersion parameters following the conventions of
*
* The variable is wavelength in micro meter
* The Lorentzians that are described by the factors is assumed to be of the form
*
*           a*x^2
* L(x) = -----------
*         x^2 - b
* where a is named preFactor in the JSON file and b is named resonance.
*/
class SellmeierMaterial
{
public:
  SellmeierMaterial(){};

  /** Loads material parameters from JSON file */
  void load( const char* fname );

  /** Computes the dispersion parameters used in MEEP. Consult: http://ab-initio.mit.edu/wiki/index.php/Dielectric_materials_in_Meep*/
  void getMEEPLorentzian( double lengthscaleInMicroMeter, unsigned int indx, double &sigma, double &resonnanceAngularFreq ) const;
  double epsInf{1.0};

  /** Returns the number of lorentzians used in the approximation */
  unsigned int nLorentzians() const{ return lorentzians.size(); };

  /** Returns epsilon value corresponding to the given wavelength (in um)*/
  double getEpsilon( double lambda ) const;

  /** Returns epsilon corresponding to the given MEEP frequency (dimless) and lengthscale (in um) */
  double getEpsilon( double lengthscale, double MEEPangFreq ) const;
protected:
  std::string fname{""};
  std::vector<Lorentzian> lorentzians;

  /** Chech that the JSON file contains the required fields */
  void checkRequiredFields( const Json::Value &root ) const;
};
#endif
