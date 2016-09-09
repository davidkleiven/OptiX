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
//#define DEBUG
//#define PRINT_BEM_MATRIX
//#define PRINT_RHS_VECTOR
#define DOUBLE_COMPARISON_ZERO 1E-5
#define PWVAC (0.5/ZVAC)

const double PI = acos(-1.0);
enum class Polarisation_t{S, P};

using namespace std;
typedef std::complex<double> cdouble;

/**
* @brief Compute the p-polarised E field based on H=(Hx,Hy,Hz)=(1.0,0.0,0.0)
* @param[in] kHat - unit vector in which the wave propagates
* @param[out] E0 - p-polarised field amplitude
*/
void getE0_p( const double kHat[3], cdouble E0[3] )
{
  // Assuming the H field is H=(Hx, Hy, Hz) = (1.0, 0.0, 0.0)
  E0[0] = 0.0;
  E0[1] = -kHat[2];
  E0[2] = kHat[1];
}

/**
* @brief Computes the poytning vector
* @param[in] EH  field components EH={Ex,Ey,Ez,Hx,Hy,Hz}
* @param[out} poynting - the resulting time averaged Poynting vector
*/
void poyntingVector(const cdouble EH[6], double poynting[3])
{
  const cdouble *E = EH;
  const cdouble *H = EH+3;
  poynting[0] = 0.5*real(E[1]*std::conj(H[2]) - E[2]*std::conj(H[1]));
  poynting[1] = 0.5*real(E[2]*std::conj(H[0]) - E[0]*std::conj(H[2]));
  poynting[2] = 0.5*real(E[0]*std::conj(H[1]) - E[1]*std::conj(H[0]));
}

/**
* @brief Computes the amplitude of a complex vector
* @param[in] vec - complex vector 
* @return Amplitude squared of vec
*/
double getAmplitude(const cdouble vec[3])
{
  return pow(std::abs(vec[0]),2) + pow(std::abs(vec[1]),2) + pow(std::abs(vec[2]),2);
}

/**
* @brief Computes the cross product between two vectors
* @param[in] vec1 - first vector
* @param[in] vec2 - second vector
* @param[out] out = vec1 x vec2
*/
void cross(const double vec1[3], const double vec2[3], double out[3])
{
  out[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  out[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
  out[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}

/**
* @brief Computes the angle of a vector with the z-axis
* @param[in] vec
* @return Angle with the z axis in degrees
*/
double angleWithZaxis(const double vec[3])
{
  double amp = sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
  return std::acos(vec[2]/amp)*180.0/PI;
}


/**
* @brief Computes the flux through plane surface
* @param[in] poynting - poytning vector of wave
* @param[in] nHat - unit normal vector of the surface
* @return FLux through surface
*/
double flux(const double poynting[3], const double nHat[3])
{
  return poynting[0]*nHat[0] + poynting[1]*nHat[1] + poynting[2]*nHat[2];
}
  

/**
* @brief Computes cosine of angle between two vectors
* @param[in] vec1 - first vector
* @param[in] vec2 - second vector
* @return Cosine of angle between vec1 and vec2
*/
double cosAlpha( const double vec1[3], const double vec2[3] )
{
  double dotProd = 0.0;
  double absv1 = 0.0;
  double absv2 = 0.0;
  for ( unsigned int i=0;i<3;i++ )
  {
    dotProd += vec1[i]*vec2[i];
    absv1 += vec1[i]*vec1[i];
    absv2 += vec2[i]*vec2[i];
  }
  absv1 = sqrt(absv1);
  absv2 = sqrt(absv2);
  return dotProd/( absv1*absv2 );
}  

/**
* @brief Check if two vectors are parallel
* @param[in] vec1 - first vector
* @param[in] vec2 - second vector
* @return true if parallell, false if not
*/
bool isParalell( const double vec1[3], const double vec2[3] )
{
  return abs( cosAlpha(vec1,vec2) - 1.0 ) < DOUBLE_COMPARISON_ZERO;
}

/**
* @brief Check if two vectors are perpendicular
* @param[in] vec1 - first vector
* @param[in] vec2 - second vector
* @return true if perpendicular, false if not
*/
bool isPerpendicular( const double vec1[3], const double vec2[3] )
{
  return abs( cosAlpha(vec1,vec2) ) < DOUBLE_COMPARISON_ZERO;
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
  E0_s[0].real(1.0);
  E0_s[0].imag(0.0);
  E0_s[1].real(0.0);
  E0_s[1].imag(1.0);
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
  const double kR[N_runs] = {5.0};
  // Assembling BEM matrix
  const double detectorPosition = 20.0;
  const double deviationMax = 1.0*detectorPosition;
  const unsigned int nDetectorPixelsInEachDirection = 200;
  double intensity[nDetectorPixelsInEachDirection][nDetectorPixelsInEachDirection];
  
  for ( unsigned int run=0;run<N_runs;run++)
  {
    std::cout << "*************************************************************\n";
    std::cout << "Run="<<run<<std::endl;
    pw.SetnHat(kHat);

    std::cout << "Assembling BEM matrix..." << std::flush;
    geo.AssembleBEMMatrix(static_cast<cdouble>(kR[run]), matrix);
    std::cout << " done\n";

    matrix->LUFactorize();
 
    std::cout << "Assembling rhs vector..." << std::flush;
    geo.AssembleRHSVector(static_cast<cdouble>(kR[run]), &pw, rhsVec);
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
    
    // Store fields and flux
    std::cout << "Evaluating fields...\n";
    double monitorPosition[3];
    monitorPosition[2] = detectorPosition;
    for ( unsigned int i=0;i<nDetectorPixelsInEachDirection;i++)
    {
      monitorPosition[0] = -deviationMax + 2.0*deviationMax*static_cast<double>(i)/static_cast<double>(nDetectorPixelsInEachDirection-1);
      for (unsigned int j=0;j<nDetectorPixelsInEachDirection;j++)
      {
        monitorPosition[1] = -deviationMax + 2.0*deviationMax*static_cast<double>(j)/static_cast<double>(nDetectorPixelsInEachDirection-1);
        cdouble EH[6];
        geo.GetFields(&pw, rhsVec, kR[run], monitorPosition, EH);
        intensity[i][j] = pow(abs(EH[0]),2) + pow(abs(EH[1]),2) + pow(abs(EH[2]),2);
      }
    }
    // Save intensity
    stringstream ss;
    ss << "data/intesity" << run << ".bin";
    ofstream datafile(ss.str().c_str(), ios::binary);
    if ( !datafile.good() )
    {
      std::cout << "Could not open data file...\n";
      return 1;
    }
    datafile.write(reinterpret_cast<char*>(&intensity[0]), pow(nDetectorPixelsInEachDirection,2)*sizeof(double));
    datafile.close();
    std::cout << "done...\n";
  }
  
  delete matrix;
  delete rhsVec;

  return 0;
}

