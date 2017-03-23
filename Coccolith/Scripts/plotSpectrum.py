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
        self.reflected = []
        self.boxScattered = []
        self.boxReference = []
        self.sourceFluxReference = []
        self.asymmetry = []
        self.crossSectionInPx = -1.0
        self.freqmin = 0.0
        self.dfreq = 0.0
        self.Nfreq = 0.0
        self.uid = 0
        self.voxelSize = -1.0 # In nano meters

    def freqmax( self ):
        return self.freqmin + self.dfreq*self.Nfreq

    def plot( self, color="black", lw=1 ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        f = np.linspace( self.freqmin, self.freqmax(), len( self.ref ) )
        wavelength = self.voxelSize/f
        ax.plot( wavelength, self.trans/self.ref, color=color, lw=lw )
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Transmittivity (\$I_e/I_0\$)")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.yaxis.set_ticks_position("left")
        ax.xaxis.set_ticks_position("bottom")
        return fig, ax

    def asymmetryFactor( self, color="black", lw=1 ):
        if ( len(self.asymmetry) == 0 ):
            raise (Exception("Assymmetry not given!"))
        f = f = np.linspace( self.freqmin, self.freqmax(), len(self.asymmetry) )
        wavelength = self.voxelSize/f
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(wavelength, self.asymmetry, color=color, lw=lw)
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Asymmetry factor")
        return fig, ax

    def plotDiff( self ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        f = np.linspace( self.freqmin, self.freqmax(), len( self.ref) )
        ax.plot( self.voxelSize/f, 1.0-self.trans/self.ref, color="black" )
        ax.set_xlabel( "Wavelength (nm)" )
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
        intensityRatio = np.abs( (self.boxScattered-self.boxFluxRef)/self.sourceFluxReference )
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

    def reflection( self, color="black", lw=1 ):
        if ( len(self.ref) == 0 ):
            raise( Exception("Reflection not initialized!") )
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        f = np.linspace( self.freqmin, self.freqmax(), len( self.ref ) )
        wavelength = self.voxelSize/f
        backScattered = np.abs( self.reflected/self.ref )
        ax.plot( wavelength, backScattered, color=color, lw=lw )
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Reflectivity (\$I_r/I_0\$)")
        return fig, ax

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

        if ( "reflPlaneScat" in hf.keys() ):
            specPlot.reflected = np.array( hf.get("reflPlaneScat") )

        if ( "boxFluxRef" in hf.keys() ):
            specPlot.boxFluxRef = np.array( hf.get("boxFluxRef") )

        if ( "Asymmetry" in hf.keys() ):
            specPlot.asymmetry = np.array( hf.get("Asymmetry") )
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
    try:
        fig4, ax4 = specPlot.reflection()
    except Exception as exc:
        print (str(exc))

    try:
        fig5, ax5 = specPlot.asymmetryFactor()
    except Exception as exc:
        print (str(exc))
    plt.show()

if __name__ == "__main__":
    main( sys.argv[1:] )
