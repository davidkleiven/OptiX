#include "waveGuide.hpp"
#include "cladding.hpp"
#include "solver1D.hpp"
#include <stdexcept>
#include <jsoncpp/json/writer.h>
#include <H5Cpp.h>
#include <hdf5_hl.h>

using namespace std;

void WaveGuide1DSimulation::setCladding( const Cladding &newCladding )
{
  cladding = &newCladding;
}

void WaveGuide1DSimulation::setSolver( Solver1D &solv )
{
  solver = &solv;
  solver->setGuide(*this);
}

void WaveGuide1DSimulation::solve()
{
  if ( cladding == NULL )
  {
    throw( runtime_error("No cladding specified"));
  }
  else if ( solver == NULL )
  {
    throw( runtime_error("No solver specified"));
  }
  solver->solve();
}

void WaveGuide1DSimulation::save( const string &fname ) const
{
  string h5fname = fname+".h5";
  string jsonfname = fname +".json";
  string h5dataname = "solution";

  hid_t file_id;
  hsize_t dims = solver->getSolution().size();
  unsigned int rank = 1;
  file_id = H5Fcreate(h5fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

  H5LTmake_dataset( file_id, h5dataname.c_str(), rank, &dims, H5T_NATIVE_DOUBLE, &solver->getSolution()[0]);

  H5Fclose(file_id);
}
