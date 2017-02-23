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

    def freqmax( self ):
        return self.freqmin + self.dfreq*self.Nfreq

    def plot( self ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        f = np.linspace( self.freqmin, self.freqmax(), len( self.data ) )
        ax.plot( f, self.data, color="black" )
        ax.set_xlabel("\$\omega (c/L)\$")
        ax.set_ylabel("Transmittivity")
        ax.spines["right"].set_visible( False )
        ax.spines["top"].set_visible( False )
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

    fig, ax = specPlot.plot()
    plt.show()

if __name__ == "__main__":
    main( sys.argv[1:] )
