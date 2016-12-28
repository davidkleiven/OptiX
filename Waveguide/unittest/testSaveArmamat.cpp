#include "paraxialSimulation.hpp"
#include <armadillo>
#include <H5Cpp.h>
#include "controlFile.hpp"
#include <string>
#include <iostream>

using namespace std;
/** Class for testing the matrix */
class SaveTest: public ParaxialSimulation
{
public:
  SaveTest(): ParaxialSimulation("SaveTest"), matrix(6, 8){};
  arma::mat matrix;
  virtual void save( ControlFile &ctl ) override;
};

/** Implementation of the member functions */
void SaveTest::save( ControlFile &ctl )
{
  string fname("data/testArmaMatrixSave.h5");
  file = new H5::H5File( fname.c_str(), H5F_ACC_TRUNC );
  saveArmaMat( matrix, "matrix" );
  clog << "File written to " << fname << endl;
}

/** Main function */
int main()
{
  SaveTest sim;

  // Fill and print matrix
  for ( unsigned int i=0;i<sim.matrix.n_rows;i++ )
  {
    for ( unsigned int j=0;j<sim.matrix.n_cols;j++ )
    {
      sim.matrix(i,j) = i*sim.matrix.n_rows + j+10;
      cout << sim.matrix(i,j) << " ";
    }
    cout << endl;
  }

  ControlFile ctl;
  sim.save( ctl );
  return 0;
}
