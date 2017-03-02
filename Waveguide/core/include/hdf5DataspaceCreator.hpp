#ifndef HDF5_DATASPACE_CREATOR_H
#define HDF5_DATASPACE_CREATOR_H
#include <stdexcept>
#include <armadillo>
#include <vector>

struct DataspaceBase
{
  int rank{1};
};

template <class arrayType>
class DataspaceCreator: public DataspaceBase
{
public:
  void setDims( const arrayType &type, hsize_t fdim[] ){ throw (std::runtime_error("Armadillo type not supported!")); };
};

template<>
class DataspaceCreator<arma::vec>: public DataspaceBase
{
public:
  void setDims( const arma::vec &vec, hsize_t fdim[] )
  {
    fdim[0] = vec.n_elem;
    rank = 1;
  }
};

template <>
class DataspaceCreator<arma::mat>: public DataspaceBase
{
public:
  void setDims( const arma::mat &mat, hsize_t fdim[] )
  {
    fdim[0] = mat.n_cols;
    fdim[1] = mat.n_rows;
    rank = 2;
  }
};

template <>
class DataspaceCreator<arma::cube>: public DataspaceBase
{
public:
  void setDims( const arma::cube &cube, hsize_t fdim[] )
  {
    fdim[0] = cube.n_rows;
    fdim[1] = cube.n_cols;
    fdim[2] = cube.n_slices;
    rank = 3;
  }
};
#endif
