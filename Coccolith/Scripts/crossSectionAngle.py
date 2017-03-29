import sys
import numpy as np
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 18
mpl.rcParams["axes.unicode_minus"] = False
from matplotlib import pyplot as plt
import json
import h5py as h5

class ScatteringPlots:
    def __init__( self, jsonInput ):
        self.good = True
        try:
            infile = open( jsonInput, 'r' )
            self.input = json.load( infile )
            infile.close()
        except:
            print ("Error when initializint the ScatteringPlots class")
            self.good = False

    def checkh5file( self, keys ):
        if ( not "boxFluxScat" in keys ):
            raise Exception("No datasets named boxFluxScat!" )
        elif ( not "boxFluxRef" in keys ):
            raise Exception("No datasets named boxFluxRef!")
        elif ( not "reflPlaneRef" in keys ):
            raise Exception("No datasets named reflPlaneRef!")
        elif ( not "vxSizeNM" in keys ):
            raise Exception("No dataset named vxSizeNM!")
        elif ( not "geometrySourceVolume" in keys ):
            raise Exception("No dataset named geometrySourceVolume!")
        elif ( not "spectrumFreqs" in keys ):
            raise Exception("No dataset named spectrumFreqs!")
        elif ( not "Asymmetry" in keys ):
            raise( Exception("No dataset named Asymmetry!"))

    def crossSection( self ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        for entry in self.input["datasets"]:
            with h5.File( entry["file"], 'r' ) as hf:
                self.checkh5file(hf.keys())
                scatFlux = np.array( hf.get("boxFluxScat") )
                refFlux = np.array( hf.get("boxFluxRef") )
                influx = np.array( hf.get("reflPlaneRef"))
                vxsize = np.array( hf.get("vxSizeNM") )[0]
                geosize = np.array( hf.get("geometrySourceVolume") )
                freqs = np.array( hf.get("spectrumFreqs") )
                # freqs[0]: Minimum frequency
                # freqs[1]: Frequency stepsize
                # freqs[2]: Number of frequencies
                f = np.linspace( freqs[0], freqs[0]+freqs[1]*freqs[2], freqs[2] )
                dx = geosize[3]-geosize[0]
                dy = geosize[4]-geosize[1]
                dz = geosize[5]-geosize[2]
                if ( entry["propagationDirection"] == "x" ):
                    crossSection = dy*dz
                elif ( entry["propagationDirection"] == "y" ):
                    crossSection = dx*dz
                elif ( entry["propagationDirection"] == "z" ):
                    crossSection = dx*dy
                else:
                    raise Exception("Propagation direction has to x, y or z!")
                scatCross = np.abs( (scatFlux-refFlux)*vxsize**2/( influx ) )
                wavelength = vxsize/f
                ax.plot( wavelength, scatCross, color=entry["color"], label="\$\SI{%d}{\degree}\$"%(entry["angle"]))
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Scattering cross section (nm\$^2\$)")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.yaxis.set_ticks_position("left")
        ax.xaxis.set_ticks_position("bottom")
        ax.legend(loc="upper center", frameon=False, labelspacing=0.05)
        return fig, ax

    def asymmetry( self ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        for entry in self.input["datasets"]:
            with h5.File(entry["file"], 'r') as hf:
                asym = np.array( hf.get("Asymmetry") )
                freqs = np.array( hf.get("spectrumFreqs") )
                f = np.linspace( freqs[0], freqs[0]+freqs[1]*freqs[2], freqs[2] )
                vxsize = np.array( hf.get("vxSizeNM") )[0]
                wavelength = vxsize/f
                ax.plot( wavelength, asym, color=entry["color"], label="\$\SI{%d}{\degree}\$"%(entry["angle"]))
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Asymmetry factor, \$g\$")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.yaxis.set_ticks_position("left")
        ax.xaxis.set_ticks_position("bottom")
        ax.legend(loc="best", frameon=False, labelspacing=0.05)
        return fig, ax


def main( argv ):
    if ( len(argv) != 1 ):
        print ("Usage: python crossSectionAngle.py inputfile.json")
    try:
        plots = ScatteringPlots(argv[0])
        plots.crossSection()
        plots.asymmetry()
        plt.show()
    except Exception as exc:
        print (str(exc))

if __name__ == "__main__":
    main( sys.argv[1:])
