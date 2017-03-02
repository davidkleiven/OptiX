#ifndef ABSORBER_BC_H
#define ABSORBER_BC_H
#include <visa/visa.hpp>
#include <complex>
class ParaxialSimulation;

typedef std::complex<double> cdouble;

enum class ApplyDim_t{ROW,COL};

template<ApplyDim_t dim>
class DimGetter;

class Absorber
{
public:
  Absorber(){};

  /** Sets the thickness in pixels */
  void setThickness( unsigned int newThick ){ thickness = newThick; };

  /** Sets the inverse damping length in inverse pixels */
  void setInverseDampingLength( double newDamping ){ inverseDampingLength = newDamping; };

  /** Applies the exponential absorber in both ends of the signal */
  template<ApplyDim_t dim>
  void apply( arma::cx_mat &field, DimGetter<dim> &getter ) const;
private:
  double inverseDampingLength{0.0};
  unsigned int thickness{0};
};

/** Getter method for rows or cols. Default returns element mat(i,fixed)*/
template<ApplyDim_t dim>
class DimGetter
{
public:
  cdouble& operator()( unsigned int i, arma::cx_mat &mat){ return mat(i,fixed); };
  unsigned int fixed{0};
  unsigned int size( arma::cx_mat &mat ){ return mat.n_rows; };
};

/** Getter method that returns element mat(fixe,i)*/
template <>
class DimGetter<ApplyDim_t::ROW>
{
public:
  cdouble& operator()( unsigned int i, arma::cx_mat &mat ){ return mat(fixed,i); };
  unsigned int fixed{0};
  unsigned int size( arma::cx_mat &mat ){ return mat.n_cols; };
};
#endif
