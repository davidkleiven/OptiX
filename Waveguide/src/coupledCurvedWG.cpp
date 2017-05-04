#include "coupledCurvedWG.hpp"
#include "curvedWaveGuide2D.hpp"
#include "curvedWGCylCrd.hpp"
#include "cladding.hpp"
#include "controlFile.hpp"
#include <vector>
#include "solver2D.hpp"
//#include "source.hpp"
#include "gaussianBeam.hpp"
#include "planeWave.hpp"

using namespace std;

CoupledCurvedWG::CoupledCurvedWG(Coordinate_t crdsyst):WaveGuideFDSimulation("TwoCoupledWG2D"), crd(crdsyst)
{
  switch ( crd )
  {
    case Coordinate_t::CARTESIAN:
      wg1 = new CurvedWaveGuideFD();
      wg2 = new CurvedWaveGuideFD();
      break;
    case Coordinate_t::CYLINDRICAL:
      wg1 = new CurvedWGCylCrd();
      wg2 = new CurvedWGCylCrd();
      break;
   }
};

 CoupledCurvedWG::~CoupledCurvedWG()
 {
   delete wg1;
   delete wg2;
 }

bool CoupledCurvedWG::isInsideGuide( double x, double z ) const
 {
   bool inWg1 = wg1->isInsideGuide( x, z );
   bool inWg2 = false;

   if ( z > startCoupler )
   {
     inWg2 = wg2->isInsideGuide( x-separation, z );
   }
   return inWg1 || inWg2;
 }

 void CoupledCurvedWG::fillInfo( Json::Value &obj ) const
 {
   Json::Value wg1Obj;
   Json::Value wg2Obj;
   wg1->fillInfo( wg1Obj );
   wg2->fillInfo( wg2Obj );
   obj["waveguide1"] = wg1Obj;
   obj["waveguide2"] = wg2Obj;
   obj["separation"] = separation;
   obj["couplerStart"] = startCoupler;
 }
