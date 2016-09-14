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
#include <cstdlib>
#include <ctime>
//#define DEBUG
//#define PRINT_VALUES_INSIDE_LOOPS
//#define PRINT_BEM_MATRIX
//#define PRINT_RHS_VECTOR
#define DOUBLE_COMPARISON_ZERO 1E-5
#define PWVAC (0.5/ZVAC)
#define UID_DIGITS 6

const double PI = acos(-1.0);
enum class Polarisation_t {CIRCULAR, LINEAR};
const unsigned int UID_MAX = pow(10, UID_DIGITS);

using namespace std;
typedef std::complex<double> cdouble;

  
unsigned int uidFromH5filename( const string& fname )
{
  auto suffix = fname.find(".h5");
  if ( ( suffix != string::npos ) && (suffix >= UID_DIGITS ) )
  {
    string uid = fname.substr(suffix-UID_DIGITS, UID_DIGITS);
    stringstream ss;
    ss << uid;
    unsigned int int_uid;
    ss >> int_uid;
    return int_uid;
  }
  return UID_MAX;
}
  
void avgPoynting( const cdouble E[3], const cdouble H[3], double poynting[3] )
{
  poynting[0] = 0.5*real(E[1]*conj(H[2]) - E[2]*conj(H[1]));
  poynting[1] = 0.5*real(E[2]*conj(H[0]) - E[0]*conj(H[2]));
  poynting[2] = 0.5*real(E[0]*conj(H[1]) - E[1]*conj(H[0]));
}
  
void saveField( HMatrix &data, unsigned int pixels, const string &fname )
{ 
    double intensity[pixels][pixels];
    #ifdef DEBUG
      clog << "Copying files to temporary array... ";
      clog << "Number of rows in data: " << data.NR << endl;
      clog << "Number of columns in data: " << data.NC << endl;
    #endif
    for ( unsigned int i=0;i<pixels;i++)
    {
      for (unsigned int j=0;j<pixels;j++)
      {
        unsigned int indx = i*pixels+j;
        cdouble E[3];
        cdouble H[3];
        #ifdef PRINT_VALUES_INSIDE_LOOPS
          clog << "Row: " << indx << endl;
        #endif
        for ( unsigned int k=0;k<3;k++)
        {
          E[k] = data.GetEntry(indx,k);
          H[k] = data.GetEntry(indx,k+3);
          #ifdef PRINT_VALUES_INSIDE_LOOPS
            clog << "Column E: " << k << " Column H: " << k+3 << " E["<<k<<"]="<<E[k] << " H[" << k << "]=" << H[k] << endl;
          #endif
        }
        double S[3];
        avgPoynting(E,H,S);
        intensity[i][j] = pow(S[0],2) + pow(S[1],2) + pow(S[2],2);
        #ifdef PRINT_VALUES_INSIDE_LOOPS
          clog << "Sx=" << S[0] << "Sy=" << S[1] << "Sz=" << S[2] << endl;
          clog << "Intensity: " << intensity[i][j] << endl;
        #endif
      }
    }
    #ifdef DEBUG
      clog << "done\n";
    #endif

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
    std::clog << "Data written to " << fname << endl;
}

int main(int argc, char **argv)
{
  srand(time(0)); 
  unsigned int uid = rand()%UID_MAX;
  #ifdef DEBUG
    clog << "Compiled with DEBUG flag...\n";
  #endif
  #ifdef PRINT_VALUES_INSIDE_LOOPS
    clog << "Compiled with PRINT_VALUES_INSIDE_LOOPS flag...\n";
  #endif

  // Check if a solution file is specified in the input arguments
  string solutionfile("");
  for ( unsigned int i=1;i<argc;i++ )
  {
    string arg(argv[i]);
    if ( arg.find("--solution=") != string::npos )
    {
      solutionfile = arg.substr(9);
    }
  }
       
  Polarisation_t pol=Polarisation_t::LINEAR;

  string geofile("sphere.scuffgeo");
  //scuff::RWGGeometry::AssignBasisFunctionsToExteriorEdges=false;
  scuff::RWGGeometry geo = scuff::RWGGeometry(geofile.c_str());
  stringstream logfilename;
  logfilename << "sphere" << uid << ".log";
  SetLogFileName(logfilename.str().c_str());
  geo.SetLogLevel(SCUFF_VERBOSELOGGING);

  std::clog << "The geometry consists of " << geo.NumRegions << " regions\n";
  
  // Allocate BEM matrix and rhs vector
  HMatrix *matrix = geo.AllocateBEMMatrix();
  HVector *rhsVec = geo.AllocateRHSVector();

  // Source definition
  double sourcePosition[3] = {0.0,0.0,-100.0};
  double kHat[3] = {0.0,0.0,1.0};

  // Wave polarised wave
  cdouble E0_s[3];
  switch ( pol )
  {
    case Polarisation_t::CIRCULAR:
      E0_s[0].real(1.0/sqrt(2.0));
      E0_s[0].imag(0.0);
      E0_s[1].real(0.0);
      E0_s[1].imag(1.0/sqrt(2.0));
      clog << "Using circular polarised incident wave...\n";
      break;
    case Polarisation_t::LINEAR:
      E0_s[0].real(1.0);
      E0_s[0].imag(0.0);
      E0_s[1].real(0.0);
      E0_s[1].imag(0.0);
      clog << "Using linear polarised incident wave...\n";
      break;
  }  
  E0_s[2].real(0.0);
  E0_s[2].imag(0.0);
  PlaneWave pw(E0_s, kHat);
   
  for ( unsigned int i=0;i<geo.NumRegions; i++ )
  {
    double omega = 1.0; // Dummy variable as the refractive index is not freq dependent here
    std::clog << "Refractive index in region " << i << ": " << geo.RegionMPs[i]->GetRefractiveIndex(omega) << std::endl;
  }
  for ( unsigned int i=0;i<geo.NumRegions; i++ )
  {
    std::clog << "Description of region " << i << " " << geo.RegionLabels[i] << std::endl;
  }
  
  #ifdef DEBUG
    unsigned int sourceNum = 0;
    for (IncField* IF=&pw; IF != NULL; IF=IF->Next)
    {
      std::clog << "Region index of source " << sourceNum++ << ": " << IF->RegionIndex << std::endl;
    }
  #endif

  const unsigned int N_runs = 1;
  const double kR[N_runs] = {5.0};

  // Assembling BEM matrix
  const double detectorPosition = 1E4;
  const double deviationMax = 1.5*detectorPosition;
  const unsigned int nDetectorPixelsInEachDirection = 100;
  HMatrix *Xpoints = new HMatrix(nDetectorPixelsInEachDirection*nDetectorPixelsInEachDirection, 3);

  // Fill evaluation points
  for ( unsigned int i=0;i<nDetectorPixelsInEachDirection;i++ )
  {
    double x = -deviationMax + 2.0*deviationMax*static_cast<double>(i)/static_cast<double>(nDetectorPixelsInEachDirection-1);
    for ( unsigned int j=0;j<nDetectorPixelsInEachDirection;j++ )
    {
      double y = -deviationMax + 2.0*deviationMax*static_cast<double>(j)/static_cast<double>(nDetectorPixelsInEachDirection-1);
      Xpoints->SetEntry(i*nDetectorPixelsInEachDirection+j, 0, x);
      Xpoints->SetEntry(i*nDetectorPixelsInEachDirection+j, 1, y); 
      Xpoints->SetEntry(i*nDetectorPixelsInEachDirection+j, 2, detectorPosition); 
    }
  }
  
  HMatrix *evaluatedFields = new HMatrix( nDetectorPixelsInEachDirection*nDetectorPixelsInEachDirection, 6, LHM_COMPLEX );
  
  for ( unsigned int run=0;run<N_runs;run++)
  {
    uid = rand()%UID_MAX;
    double omega = kR[run];
    std::clog << "*************************************************************\n";
    std::clog << "Run="<<run<<std::endl;
    clog << "kR=" << kR[run] << endl;
    pw.SetnHat(kHat);

    if ( solutionfile == "" )
    {
      // Run new simulation
      std::clog << "Assembling BEM matrix...";
      geo.AssembleBEMMatrix(static_cast<cdouble>(omega), matrix);
      std::clog << " done\n";

      std::clog << "Computing LU decomposition... ";
      matrix->LUFactorize();
      clog << "done\n";
   
      std::clog << "Assembling rhs vector...";
      geo.AssembleRHSVector(static_cast<cdouble>(omega), &pw, rhsVec);
      std::clog << " done\n";

      #ifdef PRINT_RHS_VECTOR
        std::clog << "RHS Vector before solving...\n";
        for ( unsigned int i=0;i<rhsVec->N; i++ )
        {
          std::cout << rhsVec->GetEntry(i) << " ";
        }
        std::cout << std::endl;
      #endif

      std::clog << "Solving system of equations... " << std::flush;
      int info = matrix->LUSolve(rhsVec);
      std::clog << " done\n";
      stringstream solFname;
      solFname << "data/solution" <<uid<< ".h5";
      clog << "Exporting solution...\n";
      rhsVec->ExportToHDF5( solFname.str().c_str(), "SolutionVec" );
      clog << " done. Written to " << solFname.str() << endl;
    }
    else
    {
      // Only evaluate fields using the input file
      rhsVec->ImportFromHDF5(solutionfile.c_str(), "solution");
      uid = uidFromH5filename( solutionfile );
    }
    delete matrix;

    #ifdef PRINT_RHS_VECTOR
      std::clog << "RHS Vector after solving...\n";
      for ( unsigned int i=0;i<rhsVec->N;i++ )
      {
        std::cout << rhsVec->GetEntry(i) << " ";
      }
      std::cout << std::endl;
    #endif

    
    // Overview file
    Json::Value base;
    stringstream surfFname;
    surfFname << "data/surfaceCurrent" << uid << ".pp";
    clog << "Exporting surface currents... ";
    geo.PlotSurfaceCurrents(rhsVec, static_cast<cdouble>(omega), surfFname.str().c_str() );
    clog << " done\n";
    base["SurfaceCurrents"] = surfFname.str();

    // Store fields and flux
    std::clog << "Evaluating fields...\n";
    evaluatedFields = geo.GetFields( NULL, rhsVec, omega, Xpoints, evaluatedFields );
    stringstream ss;
    ss << "data/scattered" << uid << ".bin";

    try
    {
      saveField( *evaluatedFields, nDetectorPixelsInEachDirection, ss.str() );
      base["ScatteredField"] = ss.str();
      ss.clear();
      ss.str("");
      evaluatedFields = geo.GetFields( &pw, rhsVec, omega, Xpoints, evaluatedFields );

      ss << "data/totalfield" << uid << ".bin";
      saveField( *evaluatedFields, nDetectorPixelsInEachDirection, ss.str() );
      base["TotalField"] = ss.str();

      // Evaluate fields at a line through the center in addition
      delete Xpoints;
      delete evaluatedFields;
      Xpoints = new HMatrix(nDetectorPixelsInEachDirection, 3);
      evaluatedFields = new HMatrix(nDetectorPixelsInEachDirection, 6, LHM_COMPLEX );
      
      for ( unsigned int i=0;i<nDetectorPixelsInEachDirection;i++)
      {
        double y = -deviationMax + 2.0*deviationMax*static_cast<double>(i)/static_cast<double>(nDetectorPixelsInEachDirection-1);
        Xpoints->SetEntry(i, 0, 0.0);
        Xpoints->SetEntry(i, 1, y);
        Xpoints->SetEntry(i, 2, detectorPosition);
      }
      evaluatedFields = geo.GetFields( NULL, rhsVec, omega, Xpoints, evaluatedFields ); 
      ss.clear();
      ss.str("");
      ss << "data/xPointCenter" << uid << ".h5";
      base["XpointsCenter"] = ss.str();
      Xpoints->ExportToHDF5( ss.str().c_str(), "rVec" );
      ss.clear();
      ss.str("");
      ss << "data/fieldsCenter" << uid << ".h5";
      base["FieldCenter"] = ss.str();
      evaluatedFields->ExportToHDF5( ss.str().c_str(), "Fields" );
    }
    catch ( exception &exc )
    {
      clog << exc.what() << endl;
      return 1;
    }
    clog << "Field evaluation finished\n";
    
    cdouble eps = geo.RegionMPs[1]->GetEps(omega);
    base["Detector"]["z"] = detectorPosition;
    base["Detector"]["min"] = -deviationMax;
    base["Detector"]["max"] = deviationMax;
    base["Detector"]["pixels"] = nDetectorPixelsInEachDirection;
    base["kR"] = kR[run];
    base["eps"]["real"] = real(eps);
    base["eps"]["imag"] = imag(eps);
    base["nodes"] = rhsVec->N;
    Json::StyledWriter sw;
    ss.clear();
    ss.str("");
    ss << "data/overview" << uid << ".json";
    ofstream overview(ss.str().c_str());
    if ( !overview.good() )
    {
      clog << "Could not open file " << ss.str() << endl;
      return 1;
    }
    overview << sw.write(base) << endl;
    overview.close();
    clog << "Overview file written to " << ss.str() << endl;
  }
  
  delete rhsVec;
  delete evaluatedFields;
  delete Xpoints;

  return 0;
}

