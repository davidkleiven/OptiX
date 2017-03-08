import sys
import numpy as np
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 18
mpl.rcParams["axes.unicode_minus"] = False
from matplotlib import pyplot as plt

def main( argv ):
    if ( len(argv) != 1 ):
        print ("Usage: python3 plotDispersion.py  <datafile.csv>" )
        return 1

    data = np.loadtxt( argv[0], delimiter=",")
    wavelength = data[:,0]
    refractiveIndex = data[:,1]

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot( wavelength*1000.0, refractiveIndex, color="black" )
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Refractive index" )
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    plt.show()

if __name__ == "__main__":
    main( sys.argv[1:] )
