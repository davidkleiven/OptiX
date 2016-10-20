#include <iostream>
#include "cladding.hpp"
#include "curvedWaveGuide2D.hpp"
#include "crankNicholson.hpp"
#include "controlFile.hpp"
#include "straightWG2D.hpp"
#include "waveGuide.hpp"
#include "waveGuideRadiusCurvature.hpp"
#include "solver1D.hpp"
#include <complex>
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>
#include "harmonicInversion.hpp"
#include <jsoncpp/json/writer.h>
#include <H5Cpp.h>
#include <hdf5_hl.h>
#include <armadillo>

using namespace std;
typedef complex<double> cdouble;

bool compareCtlFiles( ControlFile &ctlEig, ControlFile &ctlFD ); // Compares that the controlfiles comes from the same system
int main( int argc, char** argv )
{
  /************************ PARSE COMMANDLINE ARGUMENTS ***********************/
  string HELP_MSG("Usage: ./eigenProject.out --ctlfd=<control file fd simulation> --ctleig=<contol file eigenmodes> [--help]");
  HELP_MSG += "help: Print this message\n";
  HELP_MSG += "ctlfd: Control file from the Finite Difference simulation (extension .json)\n";
  HELP_MSG += "ctleig: Contol file from the eigenmode simulation (extension .json)\n";
  string fnameFdCtl("");
  string fnameEigCtl("");
  for ( unsigned int i=1;i<argc;i++ )
  {
    string arg(argv[i]);
    if ( arg.find("--ctlfd=") != string::npos )
    {
      fnameFdCtl = arg.substr(8);
    }
    else if ( arg.find("--ctleig=") != string::npos )
    {
      fnameEigCtl = arg.substr(9);
    }
    else if ( arg.find("--help") != string::npos )
    {
      cout << HELP_MSG << endl;
      return 0;
    }
    else
    {
      cout << "Unknown argument " << arg << endl;
      return 0;
    }
  }

  if ( fnameFdCtl == "" )
  {
    cout << "No control file for the Finite Difference simulation specified\n";
    return 0;
  }
  else if ( fnameEigCtl == "" )
  {
    cout << "No control file for the eigenmode simulation specified\n";
    return 0;
  }
  /************** FINISHED PARSING COMMAND LINE ARGUMENTS *********************/

  // Parameters used for harmonic inversion
  unsigned int nFreq = 10;
  double freqMin = 1E-6;
  double freqMax = 0.005;
  //

  ControlFile fdctl;
  ControlFile eigctl;
  fdctl.load( fnameFdCtl );
  eigctl.load( fnameEigCtl );

  // Verify that the the parameters used are equal
  if ( !compareCtlFiles( eigctl, fdctl ) )
  {
    return 1;
  }

  CurvedWaveGuideFD wg;
  WaveGuideLargeCurvature wg1D;

  try
  {
    clog << "Loading results...\n";
    wg.init( fdctl );
    wg1D.load( eigctl );
    clog << "done\n";

    double zmax = wg.longitudinalDiscretization().max;
    double stepZ = wg.longitudinalDiscretization().step;
    double zmin = wg.longitudinalDiscretization().min;
    unsigned int Nz = (zmax-zmin)/stepZ - 1;
    arma::mat coeff(Nz, wg1D.getSolver()->getNmodes());
    for ( unsigned int mode=0; mode < wg1D.getSolver()->getNmodes(); mode++ )
    {
      cout << mode << ": " << wg1D.getSolver()->normEigenvec( mode ) << endl;
      double norm = wg1D.getSolver()->normEigenvec( mode );
      clog << "Solving for coefficient " << mode+1 << " of " << wg1D.getSolver()->getNmodes() << endl;
      for ( unsigned int iz=0; iz<Nz; iz++)
      {
        double z = zmin + iz*stepZ;
        coeff(iz,mode) = wg.project( z, wg1D, mode )/norm;
      }
    }

    stringstream fname;
    fname << "data/projectionCoeff_FDUID" << fdctl.get()["UID"].asInt() << "_eigUID" << eigctl.get()["UID"].asInt() << ".h5";
    hid_t file_id = H5Fcreate(fname.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    hsize_t dims[2] = {coeff.n_cols, coeff.n_rows};
    H5LTmake_dataset( file_id, "coefficients", 2, dims, H5T_NATIVE_DOUBLE, coeff.memptr());
    int fduid = fdctl.get()["UID"].asInt();
    int eiguid = eigctl.get()["UID"].asInt();
    H5LTset_attribute_int( file_id, "coefficients", "FDuid", &fduid, 1);
    H5LTset_attribute_int( file_id, "coefficients", "EIGuid", &eiguid, 1);
    H5LTset_attribute_double( file_id, "coefficients", "z0", &zmin, 1);
    H5LTset_attribute_double( file_id, "coefficients", "z1", &zmax, 1);
    H5Fclose(file_id);
    clog << "Projection coefficients written to " << fname.str() << endl;

    // Implement the harmonic inversion of the resulting signal
    HarmonicInversion harm;
    harm.setFreq( freqMin, freqMax, nFreq );
    for ( unsigned int mode=0;mode<wg1D.getSolver()->getNmodes(); mode++ )
    {
      vector<cdouble> complxData;
      for ( unsigned int i=0;i<coeff.n_rows; i++ )
      {
        complxData.push_back( coeff(i, mode) );
      }

      harm.solve( complxData );
      harm.addAttribute( "spacing", stepZ, mode );
    }
    stringstream harmfname;
    harmfname << "data/projectionHarminv_FDUID" << fdctl.get()["UID"].asInt() << "_eigUID" << eigctl.get()["UID"].asInt() << ".h5";
    harm.save( harmfname.str().c_str() );
  }
  catch ( exception &exc )
  {
    cout << exc.what() << endl;
    return 1;
  }
  catch (...)
  {
    cout << "An un recognized exception occured...\n";
    return 1;
  }
  return 0;
}

bool compareCtlFiles( ControlFile &ctlEig, ControlFile &ctlFD )
{
  double zero = 1E-10;
  // TODO: Make it possible to compare the delta and the beta. ---> get rid of the use of electron density and scattering length
  /*
  double deltaDiff = ctl1.get()["waveguide"]["Cladding"]["delta"].asDouble()-ctl2.get()["waveguide"]["Cladding"]["delta"].asDouble();
  if ( abs(deltaDiff) > zero )
  {
    cout << "The delta of the two files does not match. Difference: " << deltaDiff << endl;
    return false;
  }

  double betaDiff = ctl1.get()["waveguide"]["Cladding"]["beta"].asDouble()-ctl2.get()["waveguide"]["Cladding"]["beta"].asDouble();
  if ( abs(betaDiff) > zero )
  {
    cout << "The beta of the two files does not match. Difference: " << betaDiff << endl;
    return false;
  }
  */

  double widthDiff = ctlEig.get()["width"].asDouble() - ctlFD.get()["waveguide"]["Width"].asDouble();
  if ( abs(widthDiff) > zero )
  {
    cout << "The width of the waveguides does not match. Difference: " << widthDiff << endl;
    return false;
  }

  double Rdiff = ctlEig.get()["outerRadius"].asDouble() - ctlFD.get()["waveguide"]["RadiusOfCurvature"].asDouble();
  if ( abs(Rdiff) > zero )
  {
    cout << "The radius of curvature does not match. Difference: " << Rdiff << endl;
    return false;
  }
  return true;
}
