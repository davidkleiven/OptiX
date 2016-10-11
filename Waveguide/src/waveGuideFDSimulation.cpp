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
  if ( solver == NULL )
  {
    throw ( runtime_error("No solver specified!\n"));
  }

  string h5fname = fname+".h5";
  string jsonfname = fname+".json";

  hsize_t dims[2] = {nodeNumberLongitudinal(), nodeNumberTransverse()};
  unsigned int rank = 2;
  double **solution = allocateSolutionMatrix();

  solver->realPart(solution);
  hid_t file_id = H5Fcreate(h5fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  H5LTmake_dataset( file_id, "real", rank, dims, H5T_NATIVE_DOUBLE, solution);
  solver->imagPart(solution);
  H5LTmake_dataset( file_id, "imag", rank, dims, H5T_NATIVE_DOUBLE, solution);
  H5Fclose(file_id);
  deallocateSolutionMatrix(solution);
  clog << "Solution written to " << h5fname << endl;

  Json::Value base;
  Json::Value wginfo;
  Json::Value solverInfo;
  wginfo["Cladding"]["real"] = cladding->getRefractiveIndex().real();
  wginfo["Cladding"]["imag"] = cladding->getRefractiveIndex().imag();
  base["datafile"] = h5fname;
  base["name"] = name;
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

double** WaveGuideFDSimulation::allocateSolutionMatrix() const
{
  double **solution = new double*[nodeNumberLongitudinal()];
  for ( unsigned int i=0;i<nodeNumberLongitudinal();i++ )
  {
    solution[i] = new double[nodeNumberTransverse()];
  }
  return solution;
}

void WaveGuideFDSimulation::deallocateSolutionMatrix( double **matrix ) const
{
  for ( unsigned int i=0;i<nodeNumberLongitudinal();i++ )
  {
    delete [] matrix[i];
  }
  delete [] matrix;
}
