#include <cstdlib>
#include <iostream>
#include <libscuff.h>
#include <complex>
#include <cmath>
#include <jsoncpp/json/writer.h>
#include <fstream>
#include <string>

const double PI = acos(-1.0);
enum class Polarisation_t{S, P};

void getE0_p( const double kHat[3], std::complex<double> E0[3] )
{
  // Assuming the H field is H=(Hx, Hy, Hz) = (1.0, 0.0, 0.0)
  E0[0] = 0.0;
  E0[1] = -kHat[2];
  E0[2] = kHat[1];
}

void poyntingVector(const std::complex<double> EH[6], double poynting[3])
{
  poynting[0] = 0.5*real(EH[1]*std::conj(EH[5]) - EH[2]*std::conj(EH[4]));
  poynting[1] = 0.5*real(EH[2]*std::conj(EH[3]) - EH[0]*std::conj(EH[5]));
  poynting[2] = 0.5*real(EH[0]*std::conj(EH[4]) - EH[1]*std::conj(EH[3]));
}

double flux(double poynting[3], double nHat[3])
{
  return poynting[0]*nHat[0] + poynting[1]*nHat[1] + poynting[2]*nHat[2];
}
  
double getAmplitude(const std::complex<double> vec[3])
{
  return pow(std::abs(vec[0]),2) + pow(std::abs(vec[1]),2) + pow(std::abs(vec[2]),2);
}

int main(int argc, char **argv)
{
  scuff::RWGGeometry geo = scuff::RWGGeometry("halfspace.scuffgeo");
  SetLogFileName("fresnel.log");
  geo.SetLogLevel(SCUFF_VERBOSELOGGING);
  
  // Allocate BEM matrix and rhs vector
  HMatrix *matrix = geo.AllocateBEMMatrix();
  HVector *rhsVec = geo.AllocateRHSVector();

  // Source definition
  double sourcePosition[3] = {0.0,0.0,0.3};
  double kHat[3] = {0.0,0.0,-1.0};
  std::complex<double> E0_s[3] = {1.0,0.0,0.0}; 
  PlaneWave pw(E0_s, kHat);
  double omega = 1.0;
  
  // Array for storing the fields EH={Ex,Ey,Ez,Hx,Hy,Hz}
  std::complex<double> EHSource[6];
  std::complex<double> EHMonitor[6];

  double theta = 0.0;
  const double dtheta = 5.0;
  const double thetamax = 4.0;
  Polarisation_t pol[2] = {Polarisation_t::S, Polarisation_t::P};
  double kBloch[3] = {0.0,0.0,0.0};
 
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

  double monitorPosition[3]={0.0,0.0,-sourcePosition[2]};

  // Assembling BEM matrix
  while ( theta < thetamax )
  {
    angle.append(theta);
    std::cout << "Theta="<<theta<<std::endl;
    double kz = -cos(theta*PI/180.0);
    double ky = sin(theta*PI/180.0);
    kHat[1] = ky;
    kHat[2] = kz;
    kBloch[1] = ky*omega;
    pw.SetnHat(kHat);

    std::cout << "Assembling BEM matrix..." << std::flush;
    geo.AssembleBEMMatrix(omega, kBloch, matrix);
    std::cout << " done\n";

    matrix->LUFactorize();
 
    // Solve for s polarisation
    pw.SetE0(E0_s);
    std::cout << "***s-polarisation\n";
    std::cout << "Assembling rhs vector..." << std::flush;
    geo.AssembleRHSVector(omega, kBloch, &pw, rhsVec);
    std::cout << " done\n";

    std::cout << "Solving system of equations... " << std::flush;
    matrix->LUSolve(rhsVec);
    std::cout << " done\n";

    // Store fields and flux
    geo.GetFields(0, rhsVec, omega, kBloch, sourcePosition, EHSource);
    geo.GetFields(0, rhsVec, omega, kBloch, monitorPosition, EHMonitor);

    // Incident field
    std::complex<double> EHInc[6] = {1.0,0.0,0.0,0.0,kz,-ky};
    std::complex<double> EHRef[6];
   
    // Reflected field is the difference between incident
    for ( unsigned int i=0;i<6;i++ )
    {
      EHRef[i] = EHSource[i]-EHInc[i]; 
    }

    double fluxPlaneHat[3] = {0.0,0.0,1.0};

    // Compute poynting vectors
    double poyntingInc[3];
    double poyntingRef[3];
    double poyntingTrans[3];
    poyntingVector(EHInc, poyntingInc);
    poyntingVector(EHRef, poyntingRef);
    poyntingVector(EHMonitor, poyntingTrans);

    // Store values
    fluxReflected_s.append( -flux(poyntingRef, fluxPlaneHat)/flux(poyntingInc, fluxPlaneHat) ); 
    fluxTransmitted_s.append( flux(poyntingTrans, fluxPlaneHat)/flux(poyntingInc, fluxPlaneHat) );
    reflectionAmplitude_s.append( getAmplitude(EHRef)/getAmplitude(EHInc) );
    transmissionAmplitude_s.append( getAmplitude(EHMonitor)/getAmplitude(EHInc) );
    
    // Solve for p polarisation
    std::cout << "***p-polarisation\n";
    std::complex<double> E0_p[3];
    getE0_p( kHat, E0_p );
    
    pw.SetE0(E0_p);  
    std::cout << "Assembling rhs vector..." << std::flush;
    geo.AssembleRHSVector(omega, kBloch, &pw, rhsVec);
    std::cout << " done\n";

    std::cout << "Solving system of equations... " << std::flush;
    matrix->LUSolve(rhsVec);
    std::cout << " done\n";
    
    // Store fields and flux
    geo.GetFields(0, rhsVec, omega, kBloch, sourcePosition, EHSource);
    geo.GetFields(0, rhsVec, omega, kBloch, monitorPosition, EHMonitor);

    // Incident field
    EHInc[0] = E0_p[0];
    EHInc[1] = E0_p[1];
    EHInc[2] = E0_p[2];
    EHInc[3] = 1.0;
    EHInc[4] = 0.0;
    EHInc[5] = 0.0;
   
    // Reflected field is the difference between incident
    for ( unsigned int i=0;i<6;i++ )
    {
      EHRef[i] = EHSource[i]-EHInc[i]; 
    }

    // Compute poynting vectors
    poyntingVector(EHInc, poyntingInc);
    poyntingVector(EHRef, poyntingRef);
    poyntingVector(EHMonitor, poyntingTrans);

    // Store values
    fluxReflected_p.append( -flux(poyntingRef, fluxPlaneHat)/flux(poyntingInc, fluxPlaneHat) ); 
    fluxTransmitted_p.append( flux(poyntingTrans, fluxPlaneHat)/flux(poyntingInc, fluxPlaneHat) );
    reflectionAmplitude_p.append( getAmplitude(EHRef)/getAmplitude(EHInc) );
    transmissionAmplitude_p.append( getAmplitude(EHMonitor)/getAmplitude(EHInc) );
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
  spectra["AmplitudeReflected"]["p"] = transmissionAmplitude_p;

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
