#include <cstdlib>
#include <iostream>
#include <libscuff.h>
#include <complex>
#include <cmath>
#include <jsoncpp/json/writer.h>
#include <fstream>
#include <string>
#include <cassert>
#include <set>
#include <stdexcept>
//#define DEBUG
//#define PRINT_BEM_MATRIX
//#define PRINT_RHS_VECTOR
#define DOUBLE_COMPARISON_ZERO 1E-5
#define PWVAC (0.5/ZVAC)

const double PI = acos(-1.0);
enum class Polarisation_t{S, P};

using namespace std;
typedef std::complex<double> cdouble;

  
void saveField( HMatrix &data, unsigned int pixels, const string &fname )
{ 
    double intensity[pixels][pixels];
    for ( unsigned int i=0;i<pixels;i++)
    {
      for (unsigned int j=0;j<pixels;j++)
      {
        unsigned int indx = i*pixels+j;
        intensity[i][j] = pow(abs(data.GetEntry(indx,0)),2) + pow(abs(data.GetEntry(indx,1)),2) + \
                          pow(abs(data.GetEntry(indx,2)),2);
      }
    }

    // Save intensity
    ofstream datafile(fname.c_str(), ios::binary);
    if ( !datafile.good() )
    {
      string msg("Could not open file ");
      msg += fname;
      throw( runtime_error(msg) );
    }

    datafile.write(reinterpret_cast<char*>(&intensity[0]), pow(pixels,2)*sizeof(double));
    datafile.close();
    std::cout << "done...\n";
    std::cout << "Data written to " << fname << endl;
}

int main(int argc, char **argv)
{
 
  string geofile("sphere.scuffgeo");
  //scuff::RWGGeometry::AssignBasisFunctionsToExteriorEdges=false;
  scuff::RWGGeometry geo = scuff::RWGGeometry(geofile.c_str());
  SetLogFileName("sphere.log");
  geo.SetLogLevel(SCUFF_VERBOSELOGGING);

  std::cout << "The geometry consists of " << geo.NumRegions << " regions\n";
  
  // Allocate BEM matrix and rhs vector
  HMatrix *matrix = geo.AllocateBEMMatrix();
  HVector *rhsVec = geo.AllocateRHSVector();

  // Source definition
  double sourcePosition[3] = {0.0,0.0,-10.0};
  double kHat[3] = {0.0,0.0,1.0};

  // Circular polarised wave
  cdouble E0_s[3];
  E0_s[0].real(1.0/sqrt(2.0));
  E0_s[0].imag(0.0);
  E0_s[1].real(0.0);
  E0_s[1].imag(1.0/sqrt(2.0));
  E0_s[2].real(0.0);
  E0_s[2].imag(0.0);
  PlaneWave pw(E0_s, kHat);
   
  for ( unsigned int i=0;i<geo.NumRegions; i++ )
  {
    double omega = 1.0;
    std::cout << "Refractive index in region " << i << ": " << geo.RegionMPs[i]->GetRefractiveIndex(omega) << std::endl;
  }
  for ( unsigned int i=0;i<geo.NumRegions; i++ )
  {
    std::cout << "Description of region " << i << " " << geo.RegionLabels[i] << std::endl;
  }
  
  #ifdef DEBUG
    unsigned int sourceNum = 0;
    for (IncField* IF=&pw; IF != NULL; IF=IF->Next)
    {
      std::cout << "Region index of source " << sourceNum++ << ": " << IF->RegionIndex << std::endl;
    }
  #endif

  const unsigned int N_runs = 1;
  const double kR[N_runs] = {10.0};

  // Assembling BEM matrix
  const double detectorPosition = 10000.0;
  const double deviationMax = 0.3*detectorPosition;
  const unsigned int nDetectorPixelsInEachDirection = 20;
  HMatrix Xpoints(nDetectorPixelsInEachDirection*nDetectorPixelsInEachDirection, 3);

  // Fill evaluation points
  for ( unsigned int i=0;i<nDetectorPixelsInEachDirection;i++ )
  {
    double x = -deviationMax + 2.0*deviationMax*static_cast<double>(i)/static_cast<double>(nDetectorPixelsInEachDirection-1);
    for ( unsigned int j=0;j<nDetectorPixelsInEachDirection;j++ )
    {
      double y = -deviationMax + 2.0*deviationMax*static_cast<double>(j)/static_cast<double>(nDetectorPixelsInEachDirection-1);
      Xpoints.SetEntry(i*nDetectorPixelsInEachDirection+j, 0, x);
      Xpoints.SetEntry(i*nDetectorPixelsInEachDirection+j, 1, y); 
      Xpoints.SetEntry(i*nDetectorPixelsInEachDirection+j, 2, detectorPosition); 
    }
  }
  
  HMatrix evaluatedFields( nDetectorPixelsInEachDirection*nDetectorPixelsInEachDirection, 6, LHM_COMPLEX );
  
  for ( unsigned int run=0;run<N_runs;run++)
  {
    double omega = kR[run];
    std::cout << "*************************************************************\n";
    std::cout << "Run="<<run<<std::endl;
    pw.SetnHat(kHat);

    std::cout << "Assembling BEM matrix..." << std::flush;
    geo.AssembleBEMMatrix(static_cast<cdouble>(omega), matrix);
    std::cout << " done\n";

    matrix->LUFactorize();
 
    std::cout << "Assembling rhs vector..." << std::flush;
    geo.AssembleRHSVector(static_cast<cdouble>(omega), &pw, rhsVec);
    std::cout << " done\n";

    #ifdef PRINT_RHS_VECTOR
      std::cout << "RHS Vector before solving...\n";
      for ( unsigned int i=0;i<rhsVec->N; i++ )
      {
        std::cout << rhsVec->GetEntry(i) << " ";
      }
      std::cout << std::endl;
    #endif

    std::cout << "Solving system of equations... " << std::flush;
    int info = matrix->LUSolve(rhsVec);
    std::cout << " done\n";

    #ifdef PRINT_RHS_VECTOR
      std::cout << "RHS Vector after solving...\n";
      for ( unsigned int i=0;i<rhsVec->N;i++ )
      {
        std::cout << rhsVec->GetEntry(i) << " ";
      }
      std::cout << std::endl;
    #endif
    
    // Overview file
    Json::Value base;
    string surfaceFname("data/surfaceCurrent.pp");
    cout << "Exporting surface currents... ";
    geo.PlotSurfaceCurrents(rhsVec, static_cast<cdouble>(omega), surfaceFname.c_str() );
    cout << " done\n";
    base["SurfaceCurrents"] = surfaceFname;

    // Store fields and flux
    std::cerr << "Evaluating fields... ";
    geo.GetFields( NULL, rhsVec, omega, &Xpoints, &evaluatedFields );
    stringstream ss;
    ss << "data/scattered" << run << ".bin";

    try
    {
      saveField( evaluatedFields, nDetectorPixelsInEachDirection, ss.str() );
      base["ScatteredField"] = ss.str();
      ss.clear();
      ss.str("");
      geo.GetFields( &pw, rhsVec, omega, &Xpoints, &evaluatedFields );

      ss << "data/totalfield" << run << ".bin";
      saveField( evaluatedFields, nDetectorPixelsInEachDirection, ss.str() );
      base["TotalField"] = ss.str();
    }
    catch ( exception &exc )
    {
      cout << exc.what() << endl;
      return 1;
    }
    

    base["Detector"]["z"] = detectorPosition;
    base["Detector"]["min"] = -deviationMax;
    base["Detector"]["max"] = deviationMax;
    base["Detector"]["pixels"] = nDetectorPixelsInEachDirection;
    base["kR"] = kR[run];
    Json::StyledWriter sw;
    ss.clear();
    ss.str("");
    ss << "data/overview" << run << ".json";
    ofstream overview(ss.str().c_str());
    if ( !overview.good() )
    {
      cout << "Could not open file " << ss.str() << endl;
      return 1;
    }
    overview << sw.write(base) << endl;
    overview.close();
    cout << "Overview file written to " << ss.str() << endl;
  }
  
  delete matrix;
  delete rhsVec;


  return 0;
}

