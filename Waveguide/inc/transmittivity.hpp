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
    Transmittivity( const Transmittivity &other );
    Transmittivity operator =(const Transmittivity &other );
    ~Transmittivity();

    /** Links a waveguide instance to the transmittivity computation */
    void linkWaveguide( const CurvedWaveGuideFD &wg );

    /** Compute integrated transmittiviy */
    void compute( double z );

    /** Returns the result */
    const std::vector<double>& get() const { return *transmission; };

    /** Returns the intensity at zero */
    double getIntensityAtZero() const { return intensityAtZero; };
  private:
    const CurvedWaveGuideFD* guide;
    std::vector<double> *transmission;
    double intensityAtZero{1E80};
    bool computeIntensityAtZero{true};

    /** Trapezoidal integration */
    double trapezoidalIntegrateIntensityZ( unsigned int ixStart, unsigned int ixEnd ) const;

  };
};
#endif
