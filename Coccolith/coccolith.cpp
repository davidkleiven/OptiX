#include "voxelMaterial.hpp"
#include <stdexcept>
#include <iostream>

using namespace std;

int main( int argc, char** argv )
{
  VoxelMaterial material;
  try
  {
    material.loadRaw( "data/cocco8cv4_216_182_249_253.raw" );
  }
  catch( exception &exc )
  {
    cout << exc.what() << endl;
  }
  catch( ... )
  {
    cout << "Unrecognized exception!\n";
  }
  return 0;
}
