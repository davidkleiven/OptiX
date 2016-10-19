#include "waveGuide.hpp"
#include "cladding.hpp"
#include "solver1D.hpp"
#include <stdexcept>
#include <jsoncpp/json/writer.h>
#include <H5Cpp.h>
#include <hdf5_hl.h>
#include <fstream>
#include <cmath>
#include "controlFile.hpp"

using namespace std;

const double WaveGuide1DSimulation::PI = acos(-1.0);

// Local enum for handling the loading of HDF5 files
enum class LoadStatus_t { CONTINUE, FINISHED };

void WaveGuide1DSimulation::setCladding( const Cladding &newCladding )
{
  cladding = &newCladding;
}

WaveGuide1DSimulation::~WaveGuide1DSimulation()
{
  if ( solverInitializedFromLoad )
  {
    delete solver;
  }
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

void WaveGuide1DSimulation::save( ControlFile &ctl ) const
{
  string fname = ctl.getFnameTemplate();
  string h5fname = fname+".h5";
  string potfname = fname+"_potential.h5";

  hid_t file_id;
  hsize_t dims = solver->getEigenVectorSize();
  unsigned int rank = 1;
  file_id = H5Fcreate(h5fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

  for ( unsigned int i=0;i<solver->getNmodes();i++ )
  {
    stringstream dataname;
    dataname << "mode"<<i;
    double eigval = solver->getEigenvalue(i);

    // Armadillo stores column major in buffer
    H5LTmake_dataset( file_id, dataname.str().c_str(), rank, &dims, H5T_NATIVE_DOUBLE, solver->getSolution().colptr(i));
    H5LTset_attribute_double( file_id, dataname.str().c_str(), "eigenvalue", &eigval, 1);
  }

  H5Fclose(file_id);
  double potXmin = -2.0*width;
  double potXmax = width;
  writePotentialToFile( potfname, potXmin, potXmax );

  // Write information to a json file
  Json::Value solverParameters;
  solver->fillJsonObj(solverParameters);
  ctl.get()["solver"] = solverParameters;
  ctl.get()["solutionfile"] = h5fname;
  ctl.get()["name"] = name;
  ctl.get()["potentialFname"] = potfname;
  ctl.get()["potentialLabel"] = "potential";
  ctl.get()["innerRadius"] = innerRadius;
  ctl.get()["outerRadius"] = outerRadius;
  ctl.get()["width"] = width;
  ctl.get()["wavenumber"] = wavenumber;
  ctl.get()["potentialXmin"] = potXmin;
  ctl.get()["potentialXmax"] = potXmax;

  cout << "Data written to " << h5fname << endl;
  cout << "Waveguide potential written to " << potfname << endl;
}

void WaveGuide1DSimulation::writePotentialToFile( const string &fname, double xmin, double xmax ) const
{
  unsigned int N = 100;
  vector<double> pot;
  for ( unsigned int i=0;i<N;i++ )
  {
    double x = xmin + (xmax-xmin)*static_cast<double>(i)/static_cast<double>(N-1);
    pot.push_back( potential(x) );
  }
  hid_t file_id;
  hsize_t dims = pot.size();
  unsigned int rank = 1;
  file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

  H5LTmake_dataset( file_id, "potential", rank, &dims, H5T_NATIVE_DOUBLE, &pot[0]);

  H5Fclose(file_id);
}

void WaveGuide1DSimulation::load( ControlFile &ctl )
{
    hid_t file_id = H5Fopen(ctl.get()["solutionfile"].asString().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
    if ( file_id < 0 )
    {
      throw (runtime_error("Errorn when loading HDF5 file"));
    }

    hsize_t dim;
    // Assuming the datasets are named mode1, mode2 ... and has an attribute eigenvalue
    unsigned indx = 0;
    bool dimensionsOfMatrixIsSet = false;
    unsigned int maxIter = 10000;
    for ( unsigned int iter=0;iter<maxIter;iter++ )
    {
      stringstream name;
      name << "mode" << indx;
      herr_t infostatus = H5LTget_dataset_info(file_id, name.str().c_str(), &dim, NULL, NULL );
      if ( infostatus < 0 )
      {
        string msg("Error when searching for info in dataset ");
        msg += name.str();
        throw (runtime_error(msg.c_str()));
      }

      double buffer[dim];
      herr_t readstatus = H5LTread_dataset_double(file_id, name.str().c_str(), buffer );

      if ( readstatus < 0 )
      {
        clog << "Read " << indx << " datasets from the HDF5 file\n";
        H5Fclose(file_id);
        solver->loadingFinished();
        return;
      }

      double eigval;
      herr_t attrreadstatus = H5LTget_attribute_double(file_id, name.str().c_str(), "eigenvalue", &eigval);

      if ( attrreadstatus < 0 )
      {
        string msg("Error when reading attribute from dataset ");
        msg += name.str();
        throw (runtime_error(msg.c_str()));
      }
      solver->addEigenmode( eigval, buffer, dim );
    }
    clog << "Warning: Max number of datasets reaced\n";
    H5Fclose(file_id);
}
