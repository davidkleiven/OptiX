#include <cstdlib>
#include <iostream>
#include <libscuff.h>
#include <complex>
#include <cmath>
#include <jsoncpp/json/writer.h>
#include <fstream>
#include <string>
#include <cassert>
#include "FluxIntegrator.h"
#define DEBUG
//#define PRINT_BEM_MATRIX
//#define PRINT_RHS_VECTOR
#define DOUBLE_COMPARISON_ZERO 1E-5
#define PWVAC (0.5/ZVAC)

const double PI = acos(-1.0);
enum class Polarisation_t{S, P};

typedef std::complex<double> cdouble;

void getE0_p( const double kHat[3], cdouble E0[3] )
{
  // Assuming the H field is H=(Hx, Hy, Hz) = (1.0, 0.0, 0.0)
  E0[0] = 0.0;
  E0[1] = -kHat[2];
  E0[2] = kHat[1];
}

void poyntingVector(const cdouble EH[6], double poynting[3])
{
  const cdouble *E = EH;
  const cdouble *H = EH+3;
  poynting[0] = 0.5*real(E[1]*std::conj(H[2]) - E[2]*std::conj(H[1]));
  poynting[1] = 0.5*real(E[2]*std::conj(H[0]) - E[0]*std::conj(H[2]));
  poynting[2] = 0.5*real(E[0]*std::conj(H[1]) - E[1]*std::conj(H[0]));
}

double getAmplitude(const cdouble vec[3])
{
  return pow(std::abs(vec[0]),2) + pow(std::abs(vec[1]),2) + pow(std::abs(vec[2]),2);
}

void cross(const double vec1[3], const double vec2[3], double out[3])
{
  out[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  out[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
  out[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}

double angleWithZaxis(const double vec[3])
{
  double amp = sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
  return std::acos(vec[2]/amp)*180.0/PI;
}

double transmissionAngle( const double angle, double n1, double n2 )
{
  double sinT = n1*std::sin(angle)/n2;
  if ( std::abs( sinT ) > 1.0 )
  {
    return 0.0;
  }
  return std::asin(sinT);
}

double flux(const double poynting[3], const double nHat[3])
{
  return poynting[0]*nHat[0] + poynting[1]*nHat[1] + poynting[2]*nHat[2];
}
  

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

bool isParalell( const double vec1[3], const double vec2[3] )
{
  return abs( cosAlpha(vec1,vec2) - 1.0 ) < DOUBLE_COMPARISON_ZERO;
}

bool isPerpendicular( const double vec1[3], const double vec2[3] )
{
  return abs( cosAlpha(vec1,vec2) ) < DOUBLE_COMPARISON_ZERO;
}

int main(int argc, char **argv)
{
  scuff::RWGGeometry::AssignBasisFunctionsToExteriorEdges=false;
  scuff::RWGGeometry geo = scuff::RWGGeometry("halfspace.scuffgeo");
  SetLogFileName("fresnel.log");
  geo.SetLogLevel(SCUFF_VERBOSELOGGING);

  std::cout << "The geometry consists of " << geo.NumRegions << " regions\n";
  
  // Allocate BEM matrix and rhs vector
  HMatrix *matrix = geo.AllocateBEMMatrix();
  HVector *rhsVec = geo.AllocateRHSVector();

  // Source definition
  double sourcePosition[3] = {0.5,0.5,-2.0};
  double kHat[3] = {0.0,0.0,1.0};
  cdouble E0_s[3];
  E0_s[0].real(1.0);
  E0_s[0].imag(0.0);
  E0_s[1].real(0.0);
  E0_s[1].imag(0.0);
  E0_s[2].real(0.0);
  E0_s[2].imag(0.0);
  PlaneWave pw(E0_s, kHat);
  double omega = 1.0;
  
  // Array for storing the fields EH={Ex,Ey,Ez,Hx,Hy,Hz}
  cdouble EHSource[6];
  cdouble EHMonitor[6];

  double theta = 0.0;
  const double dtheta = 40.0;
  const double thetamax = 89.0;
  Polarisation_t pol[2] = {Polarisation_t::S, Polarisation_t::P};
  double kBloch[2] = {0.0,0.0};
 
  Json::Value reflectionAmplitude_s(Json::arrayValue); 
  Json::Value reflectionPhase_s(Json::arrayValue);
  Json::Value reflectionAmplitude_p(Json::arrayValue);
  Json::Value reflectionPhase_p(Json::arrayValue);
  Json::Value transmissionAmplitude_s(Json::arrayValue);
  Json::Value transmissionPhase_s(Json::arrayValue);
  Json::Value transmissionAmplitude_p(Json::arrayValue);
  Json::Value transmissionPhase_p(Json::arrayValue);
  Json::Value fluxReflected_s(Json::arrayValue);
  Json::Value fluxReflected_p(Json::arrayValue);
  Json::Value fluxTransmitted_s(Json::arrayValue);
  Json::Value fluxTransmitted_p(Json::arrayValue);
  Json::Value angle(Json::arrayValue);

  double monitorPosition[3]={sourcePosition[0],sourcePosition[1],-sourcePosition[2]};
  std::cout << "Source point is in region " << geo.GetRegionIndex(sourcePosition) << std::endl;
  std::cout << "Monitor point is in region " << geo.GetRegionIndex(monitorPosition) << std::endl;

  for ( unsigned int i=0;i<geo.NumRegions; i++ )
  {
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

  FluxIntegrator fluxInt;
  fluxInt.setnEvalPoints( 1 );
  fluxInt.setXmin( 0.0 );
  fluxInt.setXmax( 1.0 );
  fluxInt.setYmin( 0.0 );
  fluxInt.setYmax( 1.0 );
  fluxInt.fillEvaluationXY();

  // Assembling BEM matrix
  while ( theta < thetamax )
  {
    angle.append(theta);
    std::cout << "*************************************************************\n";
    std::cout << "Theta="<<theta<<std::endl;
    double kz = cos(theta*PI/180.0);
    double ky = sin(theta*PI/180.0);
    double ksource = real(geo.RegionMPs[pw.RegionIndex]->GetRefractiveIndex(omega))*omega;
    kHat[1] = ky;
    kHat[2] = kz;
    kBloch[1] = ksource*sin(theta*PI/180.0);
    pw.SetnHat(kHat);

    std::cout << "Assembling BEM matrix..." << std::flush;
    geo.AssembleBEMMatrix(static_cast<cdouble>(omega), kBloch, matrix);
    std::cout << " done\n";

    #ifdef PRINT_BEM_MATRIX
      std::cout << "**************************************************\n";
      std::cout << "Assembled BEM matrix...\n";
      for ( unsigned int i=0;i<matrix->NC;i++ )
      {
        for ( unsigned int j=0;j<matrix->NR;j++ )
        {
          std::cout << matrix->GetEntry(i,j) << " ";
        }
        std::cout << std::endl;
      }
      std::cout << "*************************************************\n";
    #endif

    matrix->LUFactorize();
 
    // Solve for s polarisation
    pw.SetE0(E0_s);
    std::cout << "***s-polarisation\n";
    std::cout << "Assembling rhs vector..." << std::flush;
    geo.AssembleRHSVector(static_cast<cdouble>(omega), kBloch, &pw, rhsVec);
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
    cdouble EHInc[6];
    geo.GetFields(NULL, rhsVec, omega, kBloch, sourcePosition, EHSource);
    geo.GetFields(NULL, rhsVec, omega, kBloch, monitorPosition, EHMonitor);
    geo.GetFields(&pw, NULL, omega, kBloch, sourcePosition, EHInc);
    
    #ifdef DEBUG
      std::cout << "LAPACK solution info flag: " << info << std::endl;
      std::cout << "Scattered Ex: " << EHSource[0] << " " << EHSource[1] << " " << EHSource[2] << std::endl;
      std::cout << "Transmitted field: " << EHMonitor[0] << " " << EHMonitor[1] << " " << EHMonitor[2] << std::endl;
      std::cout << "Incident field: " << EHInc[0] << " " << EHInc[1] << " " << EHInc[2] << std::endl;
    #endif

  
    double fluxPlaneHat[3] = {0.0,0.0,1.0};

    // Compute poynting vectors
    double poyntingInc[3];
    double poyntingRef[3];
    double poyntingTrans[3];
    poyntingVector(EHInc, poyntingInc);
    poyntingVector(EHSource, poyntingRef);
    poyntingVector(EHMonitor, poyntingTrans);

    // Poynting checks
    assert ( isParalell( poyntingInc, kHat ) );

    #ifdef DEBUG
      std::cout << "Incident Poynting: " << poyntingInc[0] << "," << poyntingInc[1] << "," << poyntingInc[2] << std::endl;
      std::cout << "Reflected Poynting: " << poyntingRef[0] << "," << poyntingRef[1] << "," << poyntingRef[2] << std::endl;
      std::cout << "Transmitted Poynting: " << poyntingTrans[0] << "," << poyntingTrans[1] << "," << poyntingTrans[2] << std::endl;

      std::cout << "Reflected z-flux: " << flux(poyntingRef, fluxPlaneHat) << std::endl;
      std::cout << "Incident z-flux: " << flux(poyntingInc, fluxPlaneHat) << std::endl;
      std::cout << "Transmitted z-flux: " << flux(poyntingTrans, fluxPlaneHat) << std::endl;
      std::cout << "Sum: " << ( flux(poyntingTrans, fluxPlaneHat) - flux(poyntingRef, fluxPlaneHat) )/flux(poyntingInc, fluxPlaneHat);
      std::cout << std::endl;

      double kHatT[3];
      poyntingVector( EHMonitor, kHatT );
      double n1 = real(geo.RegionMPs[0]->GetRefractiveIndex(omega));
      double n2 = real(geo.RegionMPs[1]->GetRefractiveIndex(omega));
      double thetaT = angleWithZaxis( kHatT );
      double expThetaT = transmissionAngle( theta*PI/180.0, n1, n2 )*180.0/PI;
      std::cout << "Transmitted angle: " << thetaT << ". Expected: " << expThetaT << std::endl;
    #endif

    // Compute flux
    double incFlux, refFlux, transFlux;
    fluxInt.setZpos( sourcePosition[2] );
    refFlux = fluxInt.scatteredFlux( geo, *rhsVec, omega, kBloch );
    incFlux = fluxInt.incidentFlux( geo, pw, omega, kBloch );
    fluxInt.setZpos( monitorPosition[2] );
    transFlux = fluxInt.scatteredFlux( geo, *rhsVec, omega, kBloch );
    std::cout << "Test flux: " << fluxInt.incidentFlux( geo, pw, omega, kBloch );

    #ifdef DEBUG
      std::cout << "Reflected z-flux integrated: " << refFlux << std::endl;
      std::cout << "Incident z-flux integrated: " << incFlux << std::endl;
      std::cout << "Transmitted z-flux integrated: " << transFlux << std::endl;
    #endif

    // Store values
    //fluxReflected_s.append( -flux(poyntingRef, fluxPlaneHat)/flux(poyntingInc, fluxPlaneHat) ); 
    //fluxTransmitted_s.append( flux(poyntingTrans, fluxPlaneHat)/flux(poyntingInc, fluxPlaneHat) );
    fluxReflected_s.append( -refFlux/incFlux ); 
    fluxTransmitted_s.append( transFlux/incFlux );
    reflectionAmplitude_s.append( getAmplitude(EHSource)/getAmplitude(EHInc) );
    transmissionAmplitude_s.append( getAmplitude(EHMonitor)/getAmplitude(EHInc) );

    reflectionPhase_s.append( std::arg( EHSource[0]/EHInc[0] ) );
    transmissionPhase_s.append( std::arg( EHMonitor[0]/EHInc[0] ) );
    
    // Solve for p polarisation
    std::cout << "***p-polarisation\n";
    cdouble E0_p[3];
    getE0_p( kHat, E0_p );
    #ifdef DEBUG
      double E0_preal[3] = {real(E0_p[0]), real(E0_p[1]), real(E0_p[2])};
      assert ( isPerpendicular( kHat, E0_preal ) );
    #endif
    
    pw.SetE0(E0_p);  
    std::cout << "Assembling rhs vector..." << std::flush;
    geo.AssembleRHSVector(omega, kBloch, &pw, rhsVec);
    std::cout << " done\n";

    std::cout << "Solving system of equations... " << std::flush;
    matrix->LUSolve(rhsVec);
    std::cout << " done\n";
    
    // Store fields and flux
    geo.GetFields(NULL, rhsVec, omega, kBloch, sourcePosition, EHSource);
    geo.GetFields(NULL, rhsVec, omega, kBloch, monitorPosition, EHMonitor);
    geo.GetFields(&pw, NULL, omega, kBloch, sourcePosition, EHInc);

    // Compute poynting vectors
    poyntingVector(EHInc, poyntingInc);
    poyntingVector(EHSource, poyntingRef);
    poyntingVector(EHMonitor, poyntingTrans);

    // Poynting checks
    assert ( isParalell( poyntingInc, kHat ) );

    // Compute flux
    fluxInt.setZpos( sourcePosition[2] );
    incFlux = fluxInt.incidentFlux( geo, pw, omega, kBloch );
    refFlux = fluxInt.scatteredFlux( geo, *rhsVec, omega, kBloch );
    fluxInt.setZpos( monitorPosition[2] );
    transFlux = fluxInt.scatteredFlux( geo, *rhsVec, omega, kBloch );

    // Store values
    //fluxReflected_p.append( -flux(poyntingRef, fluxPlaneHat)/flux(poyntingInc, fluxPlaneHat) ); 
    //fluxTransmitted_p.append( flux(poyntingTrans, fluxPlaneHat)/flux(poyntingInc, fluxPlaneHat) );
    fluxReflected_p.append( -refFlux/incFlux ); 
    fluxTransmitted_p.append( transFlux/incFlux );
    reflectionAmplitude_p.append( getAmplitude(EHSource)/getAmplitude(EHInc) );
    transmissionAmplitude_p.append( getAmplitude(EHMonitor)/getAmplitude(EHInc) );
    reflectionPhase_p.append( std::arg( EHSource[3]/EHInc[3] ) );
    transmissionPhase_p.append( std::arg( EHMonitor[3]/EHInc[3] ) );
    theta += dtheta;
  }

  delete matrix;
  delete rhsVec;

  // Store values
  Json::Value spectra;
  spectra["IncidentAngle"] = angle;
  spectra["FluxReflected"]["s"] = fluxReflected_s;
  spectra["FluxReflected"]["p"] = fluxReflected_p;
  spectra["FluxTransmitted"]["s"] = fluxTransmitted_s;
  spectra["FluxTransmitted"]["p"] = fluxTransmitted_p;
  spectra["AmplitudeReflected"]["s"] = reflectionAmplitude_s;
  spectra["AmplitudeReflected"]["p"] = reflectionAmplitude_p;
  spectra["AmplitudeTransmitted"]["s"] = transmissionAmplitude_s;
  spectra["AmplitudeTransmitted"]["p"] = transmissionAmplitude_p;
  spectra["PhaseReflected"]["s"] = reflectionPhase_s;
  spectra["PhaseReflected"]["p"] = reflectionPhase_p;
  spectra["PhaseTransmitted"]["s"] = transmissionPhase_s;
  spectra["PhaseTransmitted"]["p"] = transmissionPhase_p;
  spectra["geometry"]["Monitor"]["x"] = monitorPosition[0];
  spectra["geometry"]["Monitor"]["y"] = monitorPosition[1];
  spectra["geometry"]["Monitor"]["z"] = monitorPosition[2];
  spectra["geometry"]["Source"]["x"] = sourcePosition[0];
  spectra["geometry"]["Source"]["y"] = sourcePosition[1];
  spectra["geometry"]["Source"]["z"] = sourcePosition[2];

  if ( geo.NumRegions >= 2 )
  {
    double epsilon1 = real(geo.RegionMPs[0]->GetEps(omega));
    double epsilon2 = real(geo.RegionMPs[1]->GetEps(omega));
    double epsilonHigh = epsilon1 > epsilon2 ? epsilon1:epsilon2;
    double epsilonLow = epsilon1 > epsilon2 ? epsilon2:epsilon1;
    spectra["geometry"]["EpsilonHigh"] = epsilonHigh;
    spectra["geometry"]["EpsilonLow"] = epsilonLow;
  }
  
  std::string resultFile("data/coefficients.json");
  Json::FastWriter fw;
  std::ofstream os(resultFile.c_str());
  if ( !os.good() )
  {
    std::cerr << "Error when opening file " << resultFile << std::endl;
    return 1;
  }
  os << fw.write( spectra ) << std::endl;
  os.close(); 
  
  return 0;
}
