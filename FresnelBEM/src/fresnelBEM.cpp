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
#include "FluxIntegrator.h"
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
* @brief Computes the teoretical transmission angle
* @param[in] angle - incident angle in radians
* @param[in] n1 - refractive index in incident medium
* @param[in] n2 - refractive index in transmission medium
* @return Transmitted angle in radians
*/
double transmissionAngle( const double angle, double n1, double n2 )
{
  double sinT = n1*std::sin(angle)/n2;
  if ( std::abs( sinT ) > 1.0 )
  {
    return 0.0;
  }
  return std::asin(sinT);
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

void visualize( double ymin, double ymax, double zmin, double zmax, unsigned int Ny, unsigned int Nz, scuff::RWGGeometry &geo, \
                double omega, double kBloch[2], IncField &IF, HVector &rhs, const char* fname )
{

  // Fill LBV array
  double LBV[2][3];
  for ( unsigned int nd=0;nd<2;nd++ )
  {
    for ( unsigned int nc=0;nc<3;nc++ )
    {
      LBV[nd][nc] = geo.LBasis->GetEntryD(nc, nd);
    }
  }
    
  // Fill evaluation points
  HMatrix Xpoints(Ny*Nz, 3);
  double z = zmin;
  double dz = (zmax-zmin)/static_cast<double>( Nz );
  double dy = (ymax-ymin)/static_cast<double>( Ny );
  for ( unsigned int iz=0;iz<Nz;iz++)
  {
    double y = ymin;
    for ( unsigned int iy=0;iy<Ny;iy++)
    {
      Xpoints.SetEntry( iz*Ny+iy, 0, 0.5);
      Xpoints.SetEntry( iz*Ny+iy, 1, y );
      Xpoints.SetEntry( iz*Ny+iy, 2, z );
      y += dy;
    }
    z += dz;
  }
  HMatrix field( Ny*Nz, 6, LHM_COMPLEX );
  
  // Get fields
  geo.GetFields( &IF, &rhs, omega, kBloch, &Xpoints, &field );

  // Copy evaluation points to Json array for easier output
  Json::Value yPoints(Json::arrayValue);
  Json::Value zPoints(Json::arrayValue);
  Json::Value Ex(Json::arrayValue);
  Json::Value Ey(Json::arrayValue);
  Json::Value Ez(Json::arrayValue);
  for ( unsigned int i=0;i<Ny*Nz;i++ )
  {
    yPoints.append( real(Xpoints.GetEntry(i,1)) );
    zPoints.append( real(Xpoints.GetEntry(i,2)) );
    Ex.append( real(field.GetEntry(i,0)) );
    Ey.append( real(field.GetEntry(i,1)) );
    Ez.append( real(field.GetEntry(i,2)) );
  }
  Json::Value base;
  base["points"]["y"] = yPoints;
  base["points"]["z"] = zPoints;
  base["field"]["x"] = Ex;
  base["field"]["y"] = Ey;
  base["field"]["z"] = Ez;
  
  Json::FastWriter fw;
  std::ofstream os( fname );
  if ( !os.good() )
  {
    std::cerr << "Error when opening file " << fname << std::endl;;
    return;
  }
  os << fw.write( base );
  os.close();
  std::cout << "Data written to file " << fname << std::endl;
}
      
      
  
int main(int argc, char **argv)
{

  // Parse input arguments
  if ( argc != 6 )
  {
    std::cout << "Usage: ./fresnelBEM.out --odir=<ddir> --theta0=<minimum angle> --theta1=<maximum angle>\n";
    std::cout << "--ntheta=<number of angles> --geofile=<geo.scuffgeo>\n";
    return 1;
  }

  // Parse arguments
  double theta = 0.0;
  double thetamax = 0.0;
  double dtheta = 5.0;
  string odir("");
  string geofile("");
  
  unsigned int ntheta = 0;
  for ( unsigned int i=1;i<argc; i++ )
  {
    string arg(argv[i]);
    stringstream ss;
    if ( arg.find( "--odir=" ) != string::npos )
    {
      odir = arg.substr(7);
    }
    else if ( arg.find("--theta0=") != string::npos )
    {
      ss.clear();
      ss << arg.substr(9);
      ss >> theta;
    }
    else if ( arg.find("--theta1=") != string::npos )
    {
      ss.clear();
      ss << arg.substr(9);
      ss >> thetamax;
    }
    else if ( arg.find("--ntheta=") != string::npos )
    {
      ss.clear();
      ss << arg.substr(9);
      ss >> ntheta;
      dtheta = (thetamax - theta)/static_cast<double>(ntheta-1);
    }
    else if ( arg.find("--geofile=") != string::npos )
    {
      geofile = arg.substr(10);
    }
    else
    {
      cout << "Unknown argument " << arg << endl;
    }
  }

  // Consitency check
  if ( odir == "" )
  {
    cout << "Did not find any out directory...\n";
    return 1;
  }
  if ( geofile == "" )
  {
    cout << "Did not fint any geometry file...\n";
  }
      
  scuff::RWGGeometry::AssignBasisFunctionsToExteriorEdges=false;
  scuff::RWGGeometry geo = scuff::RWGGeometry(geofile.c_str());
  SetLogFileName("fresnel.log");
  geo.SetLogLevel(SCUFF_VERBOSELOGGING);

  std::cout << "The geometry consists of " << geo.NumRegions << " regions\n";
  
  // Allocate BEM matrix and rhs vector
  HMatrix *matrix = geo.AllocateBEMMatrix();
  HVector *rhsVec = geo.AllocateRHSVector();

  // Source definition
  double sourcePosition[3] = {0.5,0.5,-5.0};
  double kHat[3] = {0.0,0.0,1.0};
  cdouble E0_s[3];
  E0_s[0].real(1.0);
  E0_s[0].imag(0.0);
  E0_s[1].real(0.0);
  E0_s[1].imag(0.0);
  E0_s[2].real(0.0);
  E0_s[2].imag(0.0);
  PlaneWave pw(E0_s, kHat);
  double omega = 0.01;
  
  // Array for storing the fields EH={Ex,Ey,Ez,Hx,Hy,Hz}
  cdouble EHSource[6];
  cdouble EHMonitor[6];

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
  Json::Value allAbsorption;
  Json::Value absorption(Json::arrayValue);
  Json::Value absorptionPositions(Json::arrayValue);
  Json::Value absorptionValues(Json::arrayValue);
  Json::Value absorptionEntry;

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

  // Fill absorption position
  unsigned int scatReg = geo.GetRegionIndex(monitorPosition);
  cdouble n = geo.RegionMPs[scatReg]->GetRefractiveIndex(omega);
  bool fillAbsorptionArray = false;
  if ( imag(n) > 1E-12 )
  {
    double extinction = 2.0*PI/imag(n)*omega;
    unsigned int nvalues = 20;
    double xmax = 5.0*extinction;
    double dx = xmax/static_cast<double>(nvalues);
    double x = 0.0;
    for ( unsigned int i=0;i<nvalues;i++ )
    {
      absorptionPositions.append( static_cast<double>(i)*dx );
      absorptionValues.append(0.0);
    }
    allAbsorption["position"] = absorptionPositions;
    fillAbsorptionArray = true;
  }

  FluxIntegrator fluxInt;
  fluxInt.setnEvalPoints( 1 );
  fluxInt.setXmin( 0.0 );
  fluxInt.setXmax( 1.0 );
  fluxInt.setYmin( 0.0 );
  fluxInt.setYmax( 1.0 );
  fluxInt.fillEvaluationXY();

  // Assembling BEM matrix
  for ( unsigned int angleIter=0;angleIter<ntheta;angleIter++)
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
    // Loop over polarisations
    for ( unsigned int i=0;i<2;i++ )
    {
      Polarisation_t pol;
      if ( i == 0 )
      {
        pol = Polarisation_t::S;
      }
      else
      {
        pol = Polarisation_t::P;
      } 

      cdouble E0_p[3];
      switch (pol)
      {
        case Polarisation_t::S:
          pw.SetE0(E0_s);
          std::cout << "***s-polarisation\n";
          break;
        case Polarisation_t::P:
          std::cout << "***p-polarisation\n";
          getE0_p( kHat, E0_p );
          #ifdef DEBUG
            double E0_preal[3] = {real(E0_p[0]), real(E0_p[1]), real(E0_p[2])};
            assert ( isPerpendicular( kHat, E0_preal ) );
          #endif
          pw.SetE0(E0_p);  
          break;
      }
    
        
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

      // Compute fields at distance where absorption is important
      if ( fillAbsorptionArray )
      {
        for ( unsigned int i=0;i<absorptionPositions.size();i++ )
        {
          cdouble EHfield[6];
          cdouble Einc[6];
          double evaluationPosition[3] = {0.0,0.0,absorptionPositions[i].asDouble()};
          geo.GetFields(NULL, rhsVec, omega, kBloch, evaluationPosition, EHfield);
          geo.GetFields(&pw, NULL, omega, kBloch, sourcePosition, Einc);
          double amp = getAmplitude(EHfield)/getAmplitude(Einc);
          absorptionValues[i] = amp;
        }
        absorptionEntry["angle"] = angle;
        switch (pol)
        {
          case Polarisation_t::S:
            absorptionEntry["polarisation"] = "s";
            break;
          case Polarisation_t::P:
            absorptionEntry["polarisation"] = "p";
            break;
        }
        absorptionEntry["amplitude"] = absorptionValues;
        absorption.append(absorptionEntry);
      }
          
      
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
      std::cout << "Test flux: " << fluxInt.incidentFlux( geo, pw, omega, kBloch ) << std::endl;

      #ifdef DEBUG
        std::cout << "Reflected z-flux integrated: " << refFlux << std::endl;
        std::cout << "Incident z-flux integrated: " << incFlux << std::endl;
        std::cout << "Transmitted z-flux integrated: " << transFlux << std::endl;
      #endif

      // Store values
      //fluxReflected_s.append( -flux(poyntingRef, fluxPlaneHat)/flux(poyntingInc, fluxPlaneHat) ); 
      //fluxTransmitted_s.append( flux(poyntingTrans, fluxPlaneHat)/flux(poyntingInc, fluxPlaneHat) );
      switch (pol)
      {
        case Polarisation_t::S:
          fluxReflected_s.append( -refFlux/incFlux ); 
          fluxTransmitted_s.append( transFlux/incFlux );
          reflectionAmplitude_s.append( getAmplitude(EHSource)/getAmplitude(EHInc) );
          transmissionAmplitude_s.append( getAmplitude(EHMonitor)/getAmplitude(EHInc) );

          reflectionPhase_s.append( std::arg( EHSource[0]/EHInc[0] ) );
          transmissionPhase_s.append( std::arg( EHMonitor[0]/EHInc[0] ) );
          break;
        case Polarisation_t::P:
          fluxReflected_p.append( -refFlux/incFlux ); 
          fluxTransmitted_p.append( transFlux/incFlux );
          reflectionAmplitude_p.append( getAmplitude(EHSource)/getAmplitude(EHInc) );
          transmissionAmplitude_s.append( getAmplitude(EHMonitor)/getAmplitude(EHInc) );

          reflectionPhase_p.append( std::arg( EHSource[3]/EHInc[3] ) );
          transmissionPhase_p.append( std::arg( EHMonitor[3]/EHInc[3] ) );
          break;
      }
      
      // Output files currentry for the case 40 deg only
      if ( abs( theta - 40.0 ) < DOUBLE_COMPARISON_ZERO )
      {
        double twoLambda = 2.0/ksource;
        std::cout << "Outputting field for visualisation..." << std::flush;
        string incfield = odir + "/fieldInc.json";
        string transfield = odir + "fieldTrans.json";
        visualize( 0.0, 1.0, 0.0-twoLambda, 0.0, 5, 8, geo, omega, kBloch, pw, *rhsVec, incfield.c_str());
        visualize( 0.0, 1.0, 0.0, twoLambda, 5, 8, geo, omega, kBloch, pw, *rhsVec, incfield.c_str());
        std::cout << " done\n" << std::flush;
      }
    }
    theta += dtheta;
  }

  delete matrix;
  delete rhsVec;

  // Store values
  Json::Value spectra;
  allAbsorption["absMonitor"] = absorption;
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
  spectra["absorption"] = allAbsorption;

  if ( geo.NumRegions >= 2 )
  {
    double epsilon1 = real(geo.RegionMPs[0]->GetEps(omega));
    double epsilon2 = real(geo.RegionMPs[1]->GetEps(omega));
    double epsilonHigh = epsilon1 > epsilon2 ? epsilon1:epsilon2;
    double epsilonLow = epsilon1 > epsilon2 ? epsilon2:epsilon1;
    spectra["geometry"]["EpsilonHigh"] = epsilonHigh;
    spectra["geometry"]["EpsilonLow"] = epsilonLow;
  }
  
  std::string resultFile = odir+"/coefficients.json";
  Json::FastWriter fw;
  std::ofstream os(resultFile.c_str());
  if ( !os.good() )
  {
    std::cerr << "Error when opening file " << resultFile << std::endl;
    return 1;
  }
  os << fw.write( spectra ) << std::endl;
  os.close(); 
  std::cout << "Data written to: " << resultFile << std::endl;
  
  return 0;
}
