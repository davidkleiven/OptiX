import sys
sys.path.append("../FresnelFDTD")
import mplLaTeX as ml
import matplotlib as mpl
#mpl.rcParams.update(ml.params)
CALL_OVER_SSH = True
if ( CALL_OVER_SSH ):
    mpl.use("Agg")
import numpy as np
import h5py as h5
import json
from matplotlib import pyplot as plt
from scipy import interpolate
import transmission as trans

def amplitudeView( hffile ):
    keys = hffile.keys()
    nkeys = len(keys)/4

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(0, nkeys):
        amp = np.array( hffile.get("amplitude%d"%(i)) )
        freq = np.array( hffile.get("freq%d"%(i)))
        x = np.zeros(len(amp))+i
        ax.scatter(freq, x )
    ax.set_xlabel("Wavenumber")
    ax.set_ylabel("Transverse position")
    ax.set_xscale("log")
    fname = "Figures/harminv.jpeg"
    fig.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))

def amplitudeFFT( hffile ):
    data = np.array( hffile.get( "dataset" ))

    extent = [0.0, 100.0, 0.0, 2.0*np.pi/0.5]
    plt.imshow( np.abs(data)**2, extent=extent, aspect=1.0, cmap="coolwarm")#, norm=mpl.colors.LogNorm())
    plt.xlabel("$k_z$")
    plt.ylabel("$x$ (nm)")
    plt.gca().set_aspect((extent[1]-extent[0])/(extent[3]-extent[2]))
    plt.colorbar()
    fname = "Figures/fieldFFT.jpeg"
    plt.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))

def main(argv):
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        else:
            print ("Unknown argument %s"%(arg))
            return 0

    with h5.File(fname, 'r') as hf:
        amplitudeFFT(hf)
        #amplitudeView(hf)

if __name__ == "__main__":
    main( sys.argv[1:])
