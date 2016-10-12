#include "waveGuideFDSimulation.hpp"
#include "solver2D.hpp"
#include "cladding.hpp"
#include <H5Cpp.h>
#include <hdf5_hl.h>
#include <iostream>
#include <fstream>
#include <cmath>

const double PI = acos(-1.0);

using namespace std;
WaveGuideFDSimulation::WaveGuideFDSimulation(): xDisc(new Disctretization), zDisc(new Disctretization), name(""){};

WaveGuideFDSimulation::WaveGuideFDSimulation( const char* wgname ): WaveGuideFDSimulation()
{
  name = wgname;
}

void WaveGuideFDSimulation::solve()
{
  if ( solver == NULL )
  {
    throw ( runtime_error("No solver specified!"));
  }
  solver->solve();
}

WaveGuideFDSimulation::~WaveGuideFDSimulation()
{
  delete xDisc;
  delete zDisc;
}

void WaveGuideFDSimulation::setWaveLength( double lambda )
{
  wavenumber = 2.0*PI/lambda;
}

void WaveGuideFDSimulation::setTransverseDiscretization( double xmin, double xmax, double step )
{
  xDisc->min = xmin;
  xDisc->max = xmax;
  xDisc->step = step;
}

void WaveGuideFDSimulation::setLongitudinalDiscretization( double xmin, double xmax, double step )
{
  zDisc->min = xmin;
  zDisc->max = xmax;
  zDisc->step = step;
}

unsigned int WaveGuideFDSimulation::nodeNumberTransverse() const
{
  return (xDisc->max - xDisc->min)/xDisc->step + 1.0;
}

unsigned int WaveGuideFDSimulation::nodeNumberLongitudinal() const
{
  return ( zDisc->max - zDisc->min)/zDisc->step + 1.0;
}

void WaveGuideFDSimulation::setSolver( Solver2D &solv )
{
  solver = &solv;
  solver->setGuide( *this );
}

void WaveGuideFDSimulation::setCladding( const Cladding &clad )
{
  cladding = &clad;
}

void WaveGuideFDSimulation::save( const string &fname ) const
{
  save(fname, -1.0);
}

void WaveGuideFDSimulation::save( const string &fname, double intensityThreshold ) const
{
  if ( solver == NULL )
  {
    throw ( runtime_error("No solver specified!\n"));
  }

  bool useSparse = (intensityThreshold > 0.0);
  string h5fname = fname+".h5";
  string jsonfname = fname+".json";
  string wgFname = fname+"_wg.h5";

  //arma::abs(solver->getSolution()).save(h5fname.c_str(), arma::hdf5_binary);
  if ( useSparse )
  {
    sparseSave( h5fname, intensityThreshold );
  }
  else
  {
    arma::mat absSol = arma::abs(solver->getSolution());
    absSol.save(h5fname.c_str(), arma::hdf5_binary);
  }
  clog << "Solution written to " << h5fname << endl;

  saveWG( wgFname );
  clog << "Points inside waveguide written to " << wgFname << endl;

  Json::Value base;
  Json::Value wginfo;
  Json::Value solverInfo;
  wginfo["Cladding"]["delta"] = cladding->getDelta();
  wginfo["Cladding"]["beta"] = cladding->getBeta();
  base["datafile"] = h5fname;
  base["wgfile"] = wgFname;
  base["name"] = name;
  base["sparseSave"] = useSparse;
  base["sparseThreshold"] = intensityThreshold;
  base["xDiscretization"]["min"] = xDisc->min;
  base["xDiscretization"]["max"] = xDisc->max;
  base["xDiscretization"]["step"] = xDisc->step;
  base["zDiscretization"]["min"] = zDisc->min;
  base["zDiscretization"]["max"] = zDisc->max;
  base["zDiscretization"]["step"] = zDisc->step;
  fillInfo( wginfo );
  // TODO: For some reason the next line gives a segmentation fault
  //solver->fillInfo( solverInfo );
  base["solver"] = solverInfo;
  base["waveguide"] = wginfo;

  Json::StyledWriter sw;

  ofstream out(jsonfname.c_str());
  if ( !out.good() )
  {
    cerr << "Could not open file " << jsonfname << endl;
    return;
  }

  out << sw.write(base) << endl;
  out.close();
  clog << "Information written to " << jsonfname << endl;
}

double* WaveGuideFDSimulation::allocateSolutionMatrix() const
{
  double *solution = new double[nodeNumberTransverse()*nodeNumberLongitudinal()];
  return solution;
}

void WaveGuideFDSimulation::deallocateSolutionMatrix( double *matrix ) const
{
  delete [] matrix;
}

void WaveGuideFDSimulation::saveWG( const string &fname ) const
{
  vector<double> xInside;
  vector<double> zInside;

  for ( unsigned int ix=0; ix<nodeNumberTransverse(); ix++ )
  {
    double x = xDisc->min + static_cast<double>(ix)*xDisc->step;
    for ( unsigned int iz=0; iz<nodeNumberLongitudinal(); iz++ )
    {
      double z = zDisc->min + static_cast<double>(iz)*zDisc->step;
      if ( isInsideGuide( x, z ) )
      {
        xInside.push_back(x);
        zInside.push_back(z);
      }
    }
  }

  hsize_t dim = xInside.size();
  unsigned int rank = 1;

  hid_t file_id = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  H5LTmake_dataset( file_id, "xInside", rank, &dim, H5T_NATIVE_DOUBLE, &xInside[0]);
  H5LTmake_dataset( file_id, "zInside", rank, &dim, H5T_NATIVE_DOUBLE, &zInside[0]);
  H5Fclose(file_id);
}

void WaveGuideFDSimulation::sparseSave( const string &fname, double intensityThreshold ) const
{
  vector<double> xSave;
  vector<double> zSave;
  vector<double> intensity;
  unsigned int Nx = nodeNumberTransverse();
  unsigned int Nz = nodeNumberLongitudinal();
  arma::cx_mat sol = solver->getSolution();
  for ( unsigned int ix=0;ix<Nx; ix++ )
  {
    for ( unsigned int iz=0; iz<Nz; iz++ )
    {
      if ( abs(sol(ix, iz)) > intensityThreshold )
      {
        double x = xDisc->min + static_cast<double>(ix)*xDisc->step;
        double z = zDisc->min + static_cast<double>(iz)*zDisc->step;
        xSave.push_back(x);
        zSave.push_back(z);
        intensity.push_back( abs(sol(ix,iz)));
      }
    }
  }

  hsize_t dim = xSave.size();
  hid_t file_id = H5Fcreate( fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  H5LTmake_dataset( file_id, "x", 1, &dim, H5T_NATIVE_DOUBLE, &xSave[0]);
  H5LTmake_dataset( file_id, "z", 1, &dim, H5T_NATIVE_DOUBLE, &zSave[0]);
  H5LTmake_dataset( file_id, "intensity", 1, &dim, H5T_NATIVE_DOUBLE, &intensity[0]);
  H5Fclose(file_id);
}
