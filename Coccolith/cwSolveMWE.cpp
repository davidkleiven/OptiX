#include <meep.hpp>
#include <iostream>

using namespace meep;
using namespace std;

double eps( const meep::vec&r )
{
  double radius = 1.0;
  if ( pow(r.x()-2.5,2) + pow(r.y()-2.5,2) < radius*radius )
  {
    return 2.0;
  }
  return 1.0;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv); // do this even for non-MPI Meep
  try
  {
    double resolution = 2; // pixels per distance
    grid_volume v = vol3d(5,10, 15.0, resolution); // 5x10 2d cell
    structure s(v, eps, pml(1.0));
    fields f(&s);

    double freq = 0.3, fwidth = 0.1;
    continuous_src_time src(freq, fwidth);
    f.add_point_source(Ey, src, vec(1.1, 2.3, 1.25));

    double tol = 1E-8;
    int maxiter = 10000;
    int bicgstabL = 10;
    bool success =  f.solve_cw( tol, maxiter, bicgstabL );

    if ( !success )
    {
      clog << "Solution did not converge!\n";
    }

    meep::h5file *file = f.open_h5file("data/cwMWE.h5");
    f.output_hdf5( meep::Dielectric, v.surroundings(), file, false, true );
    file->prevent_deadlock();
    f.output_hdf5(EnergyDensity, v.surroundings(), file, false, true );
    file->prevent_deadlock();
  }
  catch(...)
  {
    clog << "An exception occured...\n!";
  }
  return 0;
}
