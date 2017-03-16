import sys
sys.path.append("Scripts")
import plotSpectrum as pspec
import numpy as np
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 18
mpl.rcParams["axes.unicode_minus"]=False
from matplotlib import pyplot as plt

class Chlorophyll:
    def __init__(self):
        self.chla = []
        self.chlb = []
        self.chlc = []
        self.wavelength = []

    def load( self, fname ):
        data = np.loadtxt( fname, delimiter="," )
        self.wavelength = data[:,0]
        self.chla = data[:,1]
        self.chlb = data[:,2]
        self.chlc = data[:,3]

def main( argv ):
    if ( len(argv) != 1 ):
        print ("Usage: python crossSectionWithChlorophyll.py <hdf5 file>")
        return

    spec = pspec.initSpectrum( argv[0] )
    fig, ax = spec.scatteringCrossSection( color="#e41a1c", lw=3 )
    chloro = Chlorophyll()
    chloro.load( "Materials/chlorophyll.csv" )
    f = spec.getDimlessFreq( chloro.wavelength )
    ax2 = ax.twinx()
    ax2.plot( chloro.wavelength, chloro.chla, color="#377eb8", label="A" )
    ax2.plot( chloro.wavelength, chloro.chlb, color="#4daf4a", label="B")
    ax2.plot( chloro.wavelength, chloro.chlc, color="#984ea3", label="C")
    ax2.set_ylabel("Specific Absorption Coeff. (m\$^2\$/mg)")
    dummy = ax2.plot([],[],color="#e41a1c", lw=3, label="\$\sigma_s\$")
    ax2.legend( loc="center right", frameon=False, labelspacing=0.05 )
    plt.show()

if __name__ == "__main__":
    main( sys.argv[1:] )
