import sys
import numpy as np
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 18
mpl.rcParams["axes.unicode_minus"]=False
from matplotlib import pyplot as plt
import h5py as h5

class Spectrum:
    def __init__( self ):
        self.data = []
        self.freqmin = 0.0
        self.dfreq = 0.0
        self.Nfreq = 0.0
        self.uid = 0
        self.voxelSize = -1.0 # In nano meters

    def freqmax( self ):
        return self.freqmin + self.dfreq*self.Nfreq

    def plot( self ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        f = np.linspace( self.freqmin, self.freqmax(), len( self.data ) )
        ax.plot( f, self.data, color="black" )
        ax.set_xlabel("Frequency \$(c/L)\$")
        ax.set_ylabel("Transmittivity (\$I_e/I_0\$)")
        if ( self.voxelSize > 0.0 ):
            ax2 = ax.twiny()
            ax2.set_xlim( self.voxelSize/self.freqmin, self.voxelSize/self.freqmax() )
            ax2ticks = [self.voxelSize/ax1tick for ax1tick in ax.get_xticks()]
            ax2.set_xticks( ax2ticks )
            ax2.set_xlabel("Wavelength (nm)")
        ax.yaxis.set_ticks_position("left")
        ax.xaxis.set_ticks_position("bottom")
        return fig, ax

def main( argv ):
    if ( len(argv) != 1 ):
        print ("Usage: python plotSpectrum.py <file.h5>")
        return
    fname = argv[0]

    specPlot = Spectrum()
    with h5.File( fname, 'r' ) as hf:
        specPlot.data = np.array( hf.get("spectrumTransmitted") )/np.array( hf.get("spectrumReference"))
        freqs = np.array( hf.get("spectrumFreqs") )
        specPlot.freqmin = freqs[0]
        specPlot.dfreq = freqs[1]
        specPlot.Nfreq = freqs[2]
        if ( "voxelsizeInNanoMeter" in hf.keys() ):
            specPlot.voxelSize = np.array( hf.get("voxelsizeInNanoMeter") )[0]

    # Tweak
    specPlot.voxelSize = 21.6
    fig, ax = specPlot.plot()
    plt.show()

if __name__ == "__main__":
    main( sys.argv[1:] )
