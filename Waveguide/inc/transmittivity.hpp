#ifndef TRANSMITTIVITY_H
#define TRANSMITTIVITY_H
#include <vector>

class CurvedWaveGuideFD;
class Solver2D;

namespace post
{
  class Transmittivity
  {
  public:
    Transmittivity();
    ~Transmittivity();

    /** Links a waveguide instance to the transmittivity computation */
    void linkWaveguide( const CurvedWaveGuideFD &wg );

    /** Compute integrated transmittiviy */
    void compute( double z );

    /** Returns the result */
    const std::vector<double>& get() const { return *transmission; };
  private:
    const CurvedWaveGuideFD* guide;
    std::vector<double> *transmission;
    double intensityAtZero{1.0};
    bool computeIntensityAtZero{true};

    /** Trapezoidal integration */
    double trapezoidalIntegrateIntensityZ( unsigned int ixStart, unsigned int ixEnd ) const;

  };
};
#endif
