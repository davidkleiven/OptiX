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
        self.ref = []
        self.trans = []
        self.boxScattered = []
        self.sourceFluxReference = []
        self.crossSectionInPx = -1.0
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
        f = np.linspace( self.freqmin, self.freqmax(), len( self.ref ) )
        ax.plot( f, self.trans/self.ref, color="black" )
        ax.set_xlabel("Frequency \$(c/L)\$")
        ax.set_ylabel("Transmittivity (\$I_e/I_0\$)")
        if ( self.voxelSize > 0.0 ):
            ax2 = ax.twiny()
            #ax2.set_xlim( self.voxelSize/self.freqmin, self.voxelSize/self.freqmax() )
            ax2ticks = [self.voxelSize/ax1tick for ax1tick in ax.get_xticks()]
            ax2.set_xticks( ax2ticks )
            ax2.set_xlabel("Wavelength (nm)")
        ax.yaxis.set_ticks_position("left")
        ax.xaxis.set_ticks_position("bottom")
        return fig, ax

    def plotDiff( self ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        f = np.linspace( self.freqmin, self.freqmax(), len( self.ref) )
        ax.plot( f, 1.0-self.trans/self.ref, color="black" )
        ax.set_xlabel( "Frequency \$(c/L)\$" )
        ax.set_ylabel("\$1-T\$")
        return fig, ax

    def scatteringCrossSection( self, color="black", lw=1 ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        if ( len(self.sourceFluxReference) == 0 ):
            print ("No source flux given!")
            return fig, ax
        if ( len(self.boxScattered) == 0 ):
            print ("No scattered box flux given!")
        if ( self.crossSectionInPx < 0.0 ):
            print ("No cross section in pixels given!")

        f = np.linspace( self.freqmin, self.freqmax(), len( self.ref ) )
        intensityRatio = np.abs( self.boxScattered/self.sourceFluxReference )
        area = self.crossSectionInPx*self.voxelSize*self.voxelSize/1E6

        wavelength = self.voxelSize/f
        ax.plot( wavelength, intensityRatio*area, color=color, lw=lw)
        ax.set_ylabel("Scattering cross section (\$\SI{}{\micro\meter\squared}\$)")
        ax.set_xlabel("Wavelength (nm)")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.xaxis.set_ticks_position("bottom")
        ax.yaxis.set_ticks_position("left")
        return fig, ax

    def getDimlessFreq( self, wavelenthInnm ):
        dimlessWave = wavelenthInnm/self.voxelSize
        return 1.0/dimlessWave

def initSpectrum( fname ):
    specPlot = Spectrum()
    with h5.File( fname, 'r' ) as hf:
        specPlot.ref = np.array( hf.get("spectrumReference"))
        specPlot.trans = np.array( hf.get("spectrumTransmitted"))
        freqs = np.array( hf.get("spectrumFreqs") )
        specPlot.freqmin = freqs[0]
        specPlot.dfreq = freqs[1]
        specPlot.Nfreq = freqs[2]
        if ( "vxSizeNM" in hf.keys() ):
            specPlot.voxelSize = np.array( hf.get("vxSizeNM") )[0]

        if ( "boxFluxScat" in hf.keys() ):
            specPlot.boxScattered = np.array( hf.get("boxFluxScat") )

        if ( "reflPlaneRef" in hf.keys() ):
            specPlot.sourceFluxReference = np.array( hf.get("reflPlaneRef") )

        if ( "geometrySourceVolume" in hf.keys() ):
            val = np.array( hf.get("geometrySourceVolume") )
            dx = val[3]-val[0]
            dy = val[4]-val[1]
            dz = val[5]-val[2]
            crossSections = np.array( [dx*dy,dx*dz,dy*dz] )
            specPlot.crossSectionInPx = np.max( crossSections )
    return specPlot

def main( argv ):
    if ( len(argv) != 1 ):
        print ("Usage: python plotSpectrum.py <file.h5>")
        return
    fname = argv[0]
    specPlot = initSpectrum( fname )

    # Tweak
    specPlot.voxelSize = 21.6
    fig, ax = specPlot.plot()
    fig2, ax2 = specPlot.plotDiff()
    fig3, ax3 = specPlot.scatteringCrossSection()
    plt.show()

if __name__ == "__main__":
    main( sys.argv[1:] )
