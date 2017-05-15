#include <PaxPro/genericScattering.hpp>
#include <PaxPro/controlFile.hpp>
#include <iostream>
#include <sstream>
#include <ctime>
#include <armadillo>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>

using namespace std;

/** Class for simulating a sphere with 4 layers */
class LayeredSphere: public MaterialFunction
{
public:
  void getXrayMatProp( double x, double y, double z, double &delta, double &beta ) const override
  {
    double r = sqrt( x*x+y*y+z*z );
    if ( r < rPolymer )
    {
      delta = deltaPoly;
      beta = betaPoly;
    }
    else if ( r < rNi )
    {
      delta = deltaNi;
      beta = betaNi;
    }
    else if ( r < rAu )
    {
      delta = deltaAu;
      beta = betaAu;
    }
    else if ( r < rSiO2 )
    {
      delta = deltaSiO2;
      beta = betaSiO2;
    }
    else
    {
      delta = 0.0;
      beta = 0.0;
    }
  };

  double getRmax() const { return rSiO2; };

  double rPolymer{1541.0}; // nm
  double rNi{1597.0};      // nm
  double rAu{1630.0};      // nm
  double rSiO2{1691.0};    // nm
  double deltaPoly{5.46e-6};
  double betaPoly{1.59e-8};
  double deltaNi{3.39e-5};
  double betaNi{8.71e-7};
  double deltaAu{6.27e-5};
  double betaAu{8.04e-6};
  //double deltaSiO2{1.48e-5}; // From ETB
  //double betaSiO2{2.50e-7}; // From ETB
  double deltaSiO2{9.43E-6};
  double betaSiO2{1.60E-6};
};

class OptimizeRadii
{
public:
  OptimizeRadii( LayeredSphere &sphere, GenericScattering &simulation ):sph{&sphere}, sim(&simulation)
  {
    opt = gsl_vector_alloc(6);
    stepsize = gsl_vector_alloc(6);
    double initstep = 30.0*uniformRandomNumber();
    initstep = initstep < 0.0 ? initstep-5.0:initstep+5.0;
    gsl_vector_set_all(stepsize,initstep);

    double sigma = 10.0;
    sph->rPolymer += uniformRandomNumber()*sigma;
    sph->rNi += uniformRandomNumber()*sigma;
    sph->rAu += uniformRandomNumber()*sigma;
    sph->rSiO2 += uniformRandomNumber()*sigma;

    gsl_vector_set( opt, 0, sph->rPolymer );
    gsl_vector_set( opt, 1, sph->rNi );
    gsl_vector_set( opt, 2, sph->rAu );
    gsl_vector_set( opt, 3, sph->rSiO2 );
    gsl_vector_set( opt, 4, sph->deltaSiO2*matPropScaling ); // Use numbers in the same order of magnitude
    gsl_vector_set( opt, 5, sph->betaSiO2*matPropScaling ); // Use numbers in the same order of magnitude


    s = gsl_multimin_fminimizer_alloc (T, 6);
    minex_func.f = target;
    minex_func.params = this;
    minex_func.n = 6;
    loadExperimental( "data/experimental_scattering.csv" );

    gsl_multimin_fminimizer_set (s, &minex_func, opt, stepsize);
  }

  ~OptimizeRadii()
  {
    gsl_vector_free( opt );
    gsl_vector_free( stepsize );
    gsl_multimin_fminimizer_free (s);
  }

  double uniformRandomNumber() const
  {
    return 2.0*(static_cast<double>(rand())/RAND_MAX)-1.0;
  }

  arma::mat data;
  double qminSim{0.0};
  double qmaxSim{0.0};
  arma::vec simulated;
  unsigned int maxIter{10000};
  double upperRadius{2000.0};
  double matPropScaling{1E7};

  void loadExperimental( const char* fname )
  {
    clog << "Reading experimental data from " << fname << endl;
    data.load( fname, arma::csv_ascii );
    data.col(0) /= 1E9;
  }

  /** Target function to minimize */
  static double target( const gsl_vector *x, void *params )
  {
    // Put values in to the sphere class
    OptimizeRadii* self = static_cast<OptimizeRadii*>(params);

    self->sph->rPolymer = gsl_vector_get(x,0);
    self->sph->rNi = gsl_vector_get(x,1);
    self->sph->rAu = gsl_vector_get(x,2);
    self->sph->rSiO2 = gsl_vector_get(x,3);
    self->sph->deltaSiO2 = gsl_vector_get( x, 4 )/self->matPropScaling;
    self->sph->betaSiO2 = gsl_vector_get( x, 5 )/self->matPropScaling;

    if ( self->sph->deltaSiO2 < 0.0 ) return 1E30;
    if ( self->sph->betaSiO2 < 0.0 ) return 1E30;

    for ( unsigned int i=1;i<4;i++ )
    {
      if ( gsl_vector_get(x,i) < gsl_vector_get(x,i-1) ) return 1E30;
      if ( gsl_vector_get(x,i) > self->upperRadius ) return 1E30;
    }

    self->sim->solve();

    arma::mat tempRes;
    self->sim->getFarField( tempRes );
    self->simulated = tempRes.col(tempRes.n_cols/2);

    // Normalize
    self->simulated *= arma::sum( self->data.col(1) )/arma::sum( self->simulated );
    arma::vec qSim = arma::linspace( self->qminSim, self->qmaxSim, self->simulated.n_elem );
    arma::vec res;
    arma::interp1( qSim, self->simulated, self->data.col(0), res, "*linear", 1E-16 );
    //double chiSq =  arma::sum( arma::pow(log(res)-log(self->data.col(1)),2) );

    // Normalize
    res *= arma::sum( self->data.col(1) )/arma::sum(res);
    double chiSq =  arma::sum( arma::pow(res/self->data.col(1) - 1.0,2) );

    unsigned int filterFraction = 30;
    unsigned int filterLength = self->simulated.n_elem/filterFraction;

    /*
    arma::vec simCopy(self->simulated);
    // Filter the vector
    for ( unsigned int i=0;i<self->simulated.n_elem;i++ )
    {
      if ( i < filterLength/2 )
      {
        self->simulated(i) = arma::mean(simCopy.subvec(0,i+filterLength/2) );
      }
      else if ( i > self->simulated.n_elem-filterLength/2-1 )
      {
        self->simulated(i) = arma::mean( simCopy.subvec(i-filterLength/2,self->simulated.n_elem-1) );
      }
      else
      {
        self->simulated(i) = arma::mean( simCopy.subvec(i-filterLength/2,i+filterLength/2) );
      }
    }

    arma::interp1( qSim, self->simulated, self->data.col(0), res, "*linear", 1E-16 );
    //double chSqFiltered = arma::sum( arma::pow(log(res)-log(self->data.col(1)),2) );
    double chSqFiltered = arma::sum( arma::pow(res/self->data.col(1)-1.0,2) );
    return chSqFiltered;
    return 0.9999*chSqFiltered + 0.0001*chiSq;*/
    return chiSq;

  }

  const gsl_multimin_fminimizer_type* T{gsl_multimin_fminimizer_nmsimplex2};
  gsl_multimin_fminimizer *s{nullptr};
  gsl_vector *opt{nullptr};
  gsl_vector *stepsize{nullptr};
  gsl_multimin_function minex_func;
  LayeredSphere* sph{nullptr};
  GenericScattering *sim{nullptr};
};

//=================== MAIN FUNCTION ===========================================
int main( int argc, char **argv )
{
  srand(time(0));
  if (( argc < 4 ) || ( argc > 6 ))
  {
    cout << "Usage: ./layeredSphere.out <fft or adim> <Nx> <Nz> [--optimize --stat=statfile]\n";
    return 0;
  }
  LayeredSphere material;

  GenericScattering sim("LayeredSphere");
  sim.description = "Simulation of a sphere consisting of a layered sphere 1) PMMA 2) Ni 3) Au 4) SiO2";
  sim.supressMessages = true;
  sim.setBeamWaist( 400.0*material.getRmax() );
  sim.xmin = -1.1*material.getRmax();
  sim.xmax = 1.1*material.getRmax();
  sim.ymin = -1.1*material.getRmax();
  sim.ymax = 1.1*material.getRmax();
  sim.zmin = -1.1*material.getRmax();
  sim.zmax = 1.1*material.getRmax();
  unsigned int Nx = stoi( argv[2] );
  unsigned int Nz = stoi( argv[3] );
  sim.dx = (sim.xmax-sim.xmin)/Nx;
  sim.dy = (sim.ymax-sim.ymin)/Nx;
  sim.dz = (sim.zmax-sim.zmin)/Nz;

  // Set size of exported files
  sim.exportNx = 2048;
  sim.exportNy = 2048;

  unsigned int finalSizeX = 256;
  unsigned int finalSizeZ = 256;
  sim.downSampleX = Nx/finalSizeX;
  sim.downSampleY = sim.downSampleX;
  sim.downSampleZ = Nz/finalSizeZ;
  double maxScatAngle = 0.5;
  sim.setMaxScatteringAngle( maxScatAngle );
  sim.wavelength = 0.177; // wavelength in nm

  string solver(argv[1]);
  if ( solver == "fft" )
  {
    sim.propagator = GenericScattering::SolverType_t::FFT;
  }
  else
  {
    sim.propagator = GenericScattering::SolverType_t::ADI;
  }

  sim.FFTPadLength = 8092;
  sim.setMaterial( material );
  unsigned int UID = rand()%1000000;
  try
  {
    bool runOptimization = false;
    for ( unsigned int i=0;i<argc;i++ )
    {
      string arg(argv[i]);
      if ( arg.find("--optimize") != string::npos )
      {
        runOptimization = true;
        break;
      }
    }

    if ( runOptimization )
    {
      clog << "Trying to optimize radii to experimental dataset...\n";
      OptimizeRadii optimizer( material, sim );
      double PI = acos(-1.0);
      optimizer.qminSim = -2.0*PI*sin(0.5*maxScatAngle*PI/180.0)/sim.wavelength;
      optimizer.qmaxSim = -optimizer.qminSim;
      optimizer.matPropScaling = 1E7;
      for ( unsigned int iter=0;iter<optimizer.maxIter;iter++ )
      {
        int status = gsl_multimin_fminimizer_iterate( optimizer.s );
        if ( status ) break;
        double size = gsl_multimin_fminimizer_size (optimizer.s);
        status = gsl_multimin_test_size (size, 1e-2);

        clog << "Iter="<<iter;
        for ( unsigned int i=0;i<4;i++ )
        {
          clog << " r"<<i<<"="<<gsl_vector_get( optimizer.s->x, i );
        }
        clog << " dSiO2="<<gsl_vector_get( optimizer.s->x, 4 )/optimizer.matPropScaling << " bSiO2="<<gsl_vector_get( optimizer.s->x,5 )/optimizer.matPropScaling;

        clog << " Size="<<size << " Minimizer=" << optimizer.s->fval << endl;

        if (status == GSL_SUCCESS)
        {
          cout << "Optimization converged...\n";
          break;
        }
      }

      string ofname("");
      ofstream out;
      // Write the results to a file
      for ( unsigned int i=0;i<argc;i++ )
      {
        string arg(argv[i]);
        if ( arg.find("--stat=") != string::npos )
        {
          ofname = arg.substr(7);
        }
      }

      if ( ofname != "" )
      {
        out.open( ofname.c_str(), ios_base::app );
      }
      else
      {
        stringstream ofname;
        ofname << "data/optimimLayeredSphereParams" << UID << ".csv";
        out.open(ofname.str().c_str());
      }
      if ( !out.good() )
      {
        cout << "Could not open output file...\n";
        return 1;
      }
      out << UID << ", Rpolymer," << material.rPolymer;
      out << ",RNi," << material.rNi;
      out << ",RAu" << material.rAu;
      out << ",RSiO2" << material.rSiO2;
      out << ",Minimizer: " << optimizer.s->fval;
      out << ",dSiO2=" << material.deltaSiO2 << ",bSiO2=" << material.betaSiO2 << "\n";
      out.close();
    }
    else
    {
      sim.solve();
    }

    stringstream fname;
    fname << "data/layeredSphere" << UID << ".h5";
    // Save results
    sim.save( fname.str() );
    clog << "Filename: " << fname.str() << endl;
  }
  catch( exception &exc )
  {
    cout << exc.what() << endl;
    return 1;
  }
  catch( ... )
  {
    cout << "Unrecognized exception...\n";
    return 2;
  }
  return 0;
}
