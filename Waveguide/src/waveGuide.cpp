#include "waveGuide.hpp"
#include "cladding.hpp"
#include "solver1D.hpp"
#include <stdexcept>
#include <jsoncpp/json/writer.h>
#include <H5Cpp.h>
#include <hdf5_hl.h>
#include <fstream>
#include <cmath>

using namespace std;

const double WaveGuide1DSimulation::PI = acos(-1.0);

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

  // Write information to a json file
  Json::Value base;
  Json::Value solverParameters;
  solver->fillJsonObj(solverParameters);
  base["solver"] = solverParameters;
  base["solutionfile"] = h5fname;
  base["datasetname"] = h5dataname;

  Json::StyledWriter sw;

  ofstream out(jsonfname.c_str());

  if ( !out.good() )
  {
    cout << "Could not open " << jsonfname << endl;
    return;
  }
  out << sw.write(base);
  out.close();
  cout << "Data written to " << h5fname << endl;
  cout << "Statistics and information written to " << jsonfname << endl;
}
