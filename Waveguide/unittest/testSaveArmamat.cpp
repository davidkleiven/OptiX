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
  arma::mat matrix2;
  virtual void save( ControlFile &ctl ) override;
  void printMatrices() const;
};

/** Implementation of the member functions */
void SaveTest::save( ControlFile &ctl )
{
  matrix2 = arma::pow( matrix, 2 );
  printMatrices();
  string fname("data/testArmaMatrixSave.h5");
  file = new H5::H5File( fname.c_str(), H5F_ACC_TRUNC );
  saveArray( matrix, "matrix" );
  saveArray( matrix2, "matrix2" );
  clog << "File written to " << fname << endl;
}

void SaveTest::printMatrices() const
{
  for ( unsigned int row=0;row<matrix.n_rows;row++ )
  {
    for ( unsigned int col=0;col<matrix.n_cols;col++ )
    {
      cout << matrix(row,col) << " ";
    }
    cout << endl;
  }
  cout << endl;
  for ( unsigned int row=0;row<matrix2.n_rows;row++ )
  {
    for ( unsigned int col=0;col<matrix2.n_cols;col++ )
    {
      cout << matrix2(row,col) << " ";
    }
    cout << endl;
  }
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
      sim.matrix(i,j) = i*sim.matrix.n_cols + j+10;
    }
  }

  ControlFile ctl;
  sim.save( ctl );
  return 0;
}
