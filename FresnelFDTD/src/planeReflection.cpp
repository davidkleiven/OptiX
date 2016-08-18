#include <iostream>
#include "meep.hpp"
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>
#include <stdexcept>
#include "sincSrc.h"
#include <jsoncpp/json/writer.h>
#include "dielectricSlab.h"
#define ODIRLEN 60
//#define OUTPUT_HDF5

/**
* This program takes 3 command line arguments
* Arg 1 - out directory. Directory where the datafiles will be stored
* Arg 2 - Relative dielectriv permittivity in scattering region. The wave starts in a region with epsilon=1.
* Arg 3 - angle of incidence in degrees
* Arg 4 - Polarization (s or p)
*
* Example:
*   ./planeReflection data/planewave 2.25 30 s
*/
using namespace std;

const double PI = acos(-1.0);
//----------------- MAIN FUNCTION ----------------------------//
int main(int argc, char **argv)
{
  string id("[planeRef] ");
  cout << id << "This program simulates scattering of a plane wave onto a smooth surface\n";
  meep::initialize mpi(argc, argv);

  // Read command line arguments
  if ( argc != 8 )
  {
    cout << id << "Two usages of this file:\n";
    cout << id << "Usage: ./planeReflection.out <out directory> <epsInScattered> <incident angle> <polarization> <relBandwidth>\n";
    cout << id << "<number of frequencies> <resolution>\n";
    cout << id << "The following arguments were given:\n";
    for ( unsigned int i=0;i<argc;i++ )
    {
      string arg(argv[i]);
      cout << arg << endl;
    }
    return 1;
  }
    
  const char* OUTDIR = argv[1];
  double epshigh;
  stringstream ss;
  ss << argv[2];
  ss >> epshigh;
  ss.clear();
  ss << argv[3];
  double angle;
  ss >> angle;
  char polarization = argv[4][0];
  double relFwidth;
  ss.clear();
  ss << argv[5];
  ss >> relFwidth;

  unsigned int nfreq;
  ss.clear();
  ss << argv[6];
  ss >> nfreq;

  // Check that angle is within range
  const double maxAngle = 90.0;
  if ( angle > maxAngle )
  {
    cout << id << "The incident angle is too large\n";
    cout << id << "Maximum angle is " << maxAngle << endl;
    return 1;
  }
  else if ( angle < 0.0 )
  {
    cout << id << "Negative angle given. Has to be in range [0,MAX_ANGLE)\n";
    return 1;
  }

  // Check that polarization holds a valid value
  if ( (polarization != 's') && (polarization != 'p') )
  {
    cout << id << "Invalid polarization value. Has to be either s or p\n";
    cout << id << "Value given: " << polarization << endl;
    return 1;
  }

  double resolution = 10.0;
  ss.clear();
  ss << argv[7];
  ss >> resolution;
  ss.clear();

  double freq = 0.3;
  double fwidth = freq*relFwidth;

  // Initialize geometry
  DielectricSlab geometry(resolution);
  geometry.setEpsLower(epshigh);

  // Compute kx
  double k = 2.0*PI*freq;
  geometry.setKx( k*sin( angle*PI/180.0 ) );

  const double minHeight = 2.0*geometry.getPMLThickness() + 4.0;

  // Verify that the size of the domain is big enough (for debugging only)
  assert ( geometry.getYsize() > minHeight );
  
  // Initialize
  meep::gaussian_src_time src(freq, fwidth);
  meep::component fieldComp;
  try
  {
    geometry.addSourceVol();
    geometry.addStructure();
    geometry.addField();
    
    if ( polarization == 's' )
    {
      fieldComp = meep::Ez;
    }
    else
    {
      fieldComp = meep::Hz;
    }
    geometry.addSource( src, fieldComp );
  }
  catch ( std::logic_error &exc )
  {
    cout << exc.what() << endl;
    return 1;
  }
  catch (...)
  {
    cout << "An unexpected exception occured...";
    return 1;
  }
    
  geometry.getField().set_output_directory(OUTDIR);
  
  // Write dielectric function to file
  geometry.output_hdf5( meep::Dielectric ); 

  // Add DFT fluxplane
  const double fluxPlanePosY = 0.5*(geometry.getYcPlane() + geometry.getPMLThickness());
  const double fluxRefPlanePosY = 0.5*(geometry.getYcPlane() + geometry.getSourceY());
  meep::volume dftVol = meep::volume(meep::vec(0.0,fluxPlanePosY), meep::vec(geometry.getXsize()-1.0,fluxPlanePosY));
  meep::volume dftVolR = meep::volume(meep::vec(0.0,fluxRefPlanePosY), meep::vec(geometry.getXsize()-1.0,fluxRefPlanePosY)); 
  meep::dft_flux transFluxY = geometry.getField().add_dft_flux_plane(dftVol, freq-fwidth/2.0, freq+fwidth/2.0, nfreq);
  meep::dft_flux fluxYReflected = geometry.getField().add_dft_flux_plane(dftVolR, freq-fwidth/2.0, freq+fwidth/2.0, nfreq); 

  Json::Value fieldTransmittedReal(Json::arrayValue);
  Json::Value fieldReflectionReal(Json::arrayValue);
  Json::Value timepoints(Json::arrayValue);

  // Time required to propagate over the domain with the slowest speed
  double speed = 1.0/sqrt(geometry.getEpsLower());
  double tPropagate = geometry.getField().last_source_time() + 1.2*geometry.getYsize()/speed;

  // Main loop.
  double transYWidth = abs( dftVol.get_max_corner().x() - dftVol.get_min_corner().x() );
  double timeToRegisterFourier = nfreq/fwidth;

  double tEnd = timeToRegisterFourier > tPropagate ? timeToRegisterFourier:tPropagate;

  unsigned int nOut = 10; // Number of output files
  double dt = tEnd/nOut;
  double nextOutputTime = 0.0;
  string fluxXFname("fluxYReflected");

  string outdir(OUTDIR);
  if ( outdir.find("bkg") == string::npos ) 
  {
    // Did not find bkg in the directory name. And the run is not for just producing images for visualization
    //Reading background flux from file
    // Load and subtract off the background fields
    string loadFname = "dataPlane/"+fluxXFname+".h5";
    // Check if the file can be accessed by trying to open it
    ifstream in(loadFname.c_str(), ios::binary);
    if ( !in.good() )
    {
      cout << id << "File " << loadFname << " cannot be accessed. Run a background run first. bkg has to enter in the outdir path...\n";
      return 1;
    }
    in.close();

    geometry.getField().set_output_directory("dataPlane");
    fluxYReflected.load_hdf5(geometry.getField(), fluxXFname.c_str());  
    fluxYReflected.scale_dfts(-1.0);
  }
  geometry.getField().set_output_directory(OUTDIR); 

  // Setup PML monitors
  const unsigned int nPMLMonitors = 10;
  const double dx = geometry.getPMLThickness()/static_cast<double>(nPMLMonitors);
  Json::Value pmlMonitorPos( Json::arrayValue );
  Json::Value pmlMonitorMaxAmplitude( Json::arrayValue );

  // Fill the PML monitor positions
  for ( unsigned int i=0;i<nPMLMonitors;i++ )
  {
    pmlMonitorPos.append( geometry.getPMLThickness() - static_cast<double>(i)*dx );
    pmlMonitorMaxAmplitude.append( 0.0 );
  }

  while ( geometry.getField().time() < tEnd )
  {
    geometry.getField().step();

    // Get field amplitude
    complex<double> fieldAmpTrans = geometry.getField().get_field(fieldComp, meep::vec(geometry.getXsize()/2.0, fluxPlanePosY));
    complex<double> fieldAmpRefl = geometry.getField().get_field(fieldComp, meep::vec(geometry.getXsize()/2.0, fluxRefPlanePosY));

    fieldTransmittedReal.append( real(fieldAmpTrans) );
    fieldReflectionReal.append( real(fieldAmpRefl) );
    timepoints.append( geometry.getField().time() );
 
    #ifdef OUTPUT_HDF5
      if ( field.time() > nextOutputTime )
      {
        geometry.output_hdf5(fieldComp);
        nextOutputTime += dt;
      }
    #endif

    // Get field amplitude at PML positions
    for ( unsigned int i=0;i<nPMLMonitors;i++ )
    {
      meep::vec pmlMonitorTemporaryPosition( 0.0, pmlMonitorPos[i].asDouble() );
      complex<double> pmlAmplitude = geometry.getField().get_field( fieldComp, pmlMonitorTemporaryPosition);
      if ( real( pmlAmplitude ) > pmlMonitorMaxAmplitude[i].asDouble() )
      {
        pmlMonitorMaxAmplitude[i] = real( pmlAmplitude );
      }
    } 
  } 
  #ifdef OUTPUT_HDF5
    geometry.output_hdf5(fieldComp);
  #endif

  if ( outdir.find("bkg") != string::npos )
  {
    // This is the background run --> save the fields 
    geometry.getField().set_output_directory("dataPlane");
    fluxYReflected.save_hdf5(geometry.getField(), fluxXFname.c_str());
  }

  double dfreq = fwidth/static_cast<double>(nfreq-1);
  double currentFreq = freq-fwidth/2.0;
  double *transmittedFlux = transFluxY.flux();
  double *reflectedFlux = fluxYReflected.flux();
  unsigned int numberOfNonPropagating = 0;
  Json::Value flux;
  Json::Value transmitted(Json::arrayValue);
  Json::Value reflected(Json::arrayValue);
  Json::Value freqArray(Json::arrayValue);
  Json::Value angleArray(Json::arrayValue);
  for ( unsigned int i=0;i<nfreq;i++ )
  {
    double sinCurrentAngle = freq*sin(angle*PI/180.0)/currentFreq;
    double currentAngle = 0.0;
    if ( abs(sinCurrentAngle) > 1.0 )
    {
      numberOfNonPropagating++;
    }
    else
    {
      currentAngle = asin(sinCurrentAngle)*180.0/PI;
      freqArray.append(currentFreq);
      angleArray.append(currentAngle);
      transmitted.append( transmittedFlux[i]/transYWidth );
      reflected.append( reflectedFlux[i]/transYWidth );
    }
    currentFreq += dfreq;
  }
  delete [] transmittedFlux;
  delete [] reflectedFlux;

  flux["geometry"]["sourcePosition"] = geometry.getSourceY();
  flux["geometry"]["slabPosition"] = geometry.getYcPlane();
  flux["geometry"]["xsize"] = geometry.getXsize();
  flux["geometry"]["ysize"] = geometry.getYsize();
  flux["geometry"]["EpsilonHigh"] = geometry.getEpsLower();
  flux["geometry"]["EpsilonLow"] = geometry.getEpsUpper();
  flux["frequency"] = freqArray;
  flux["incidentAngle"] = angleArray;
  flux["reflected"] = reflected;
  flux["transmitted"] = transmitted;
  
  Json::FastWriter fwFlux;

  cout << id << numberOfNonPropagating << " of " << nfreq << " modes do not propagate\n";

  // Write transmitted flux to file
  string ddir(OUTDIR);
  string fluxOut("transmittedFlux.json");
  fluxOut = ddir+"/"+fluxOut; 
  ofstream os(fluxOut.c_str());
  if ( !os.good() )
  {
    cout << id << "Problem when opening file " << fluxOut << endl;
  }
  else
  {
    os << fwFlux.write(flux) << endl;
    os.close();
  }
  
  // Write the monitors to file
  Json::Value monitors;
  monitors["time"] = timepoints;
  monitors["geometry"]["sourcePosition"] = geometry.getSourceY();
  monitors["geometry"]["slabPosition"] = geometry.getYcPlane();
  monitors["geometry"]["xsize"] = geometry.getXsize();
  monitors["geometry"]["ysize"] = geometry.getYsize();
  // TODO: Bad name change epsilonhigh to epsilonlower and epsilonlow to epsilonupper
  monitors["geometry"]["EpsilonHigh"] = geometry.getEpsLower();
  monitors["geometry"]["EpsilonLow"] = geometry.getEpsUpper();
  monitors["reflected"]["position"] = fluxRefPlanePosY;
  monitors["reflected"]["real"] = fieldReflectionReal;
  monitors["transmitted"]["position"] = fluxPlanePosY;
  monitors["transmitted"]["real"] = fieldTransmittedReal;
  monitors["PMLMonitors"]["position"] = pmlMonitorPos;
  monitors["PMLMonitors"]["MaxAmplitude"] = pmlMonitorMaxAmplitude;
  if ( fieldComp == meep::Ez )
  {
    monitors["FieldComponent"] = "Ez";
  }
  else if ( fieldComp == meep::Hz )
  {
    monitors["FieldComponent"] = "Hz";
  }
  else
  {
    monitors["FieldComponent"] = "Unknown";
  }


  // Write monitor to file
  string monitorOut("realField.json");
  monitorOut = ddir + "/" + monitorOut;
  os.open(monitorOut.c_str());
  if ( !os.good() )
  {
    cerr << id << "Error when opening file " << monitorOut << endl;
    return 1;
  }

  Json::FastWriter fw;
  os << fw.write(monitors) << endl;
  os.close();
  return 0;
}
