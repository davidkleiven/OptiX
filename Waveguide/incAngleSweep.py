import sys
sys.path.append("../FresnelFDTD")
sys.path.append("../")
import mplLaTeX as ml
import matplotlib as mpl
mpl.rcParams.update(ml.params)
import numpy as np
import json
import h5py as h5
from matplotlib import pyplot as plt
import colorScheme as cs
try:
    import colormaps as cmaps
    cmap = cmaps.viridis
except:
    cmap = "viridis"

DELTA = 4.14E-5 # Salditt et al
BETA = 3.45E-6 # Salditt et al

def main( argv ):
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        elif ( arg.find("--help") != -1 ):
            print ("Usage: python incAngleSweep.py --file=<hdf5 datafile>")
            print ("file: HDF5 file created by the incidentAngleSweep application")
            print ("help: Print this message")
            return
        else:
            print ("Unknown argument %s"%(arg))
            return

    with h5.File( fname, 'r' ) as hf:
        dset = hf.get("intensity")
        thetaMin = dset.attrs.get("thetaMin")
        thetaMax = dset.attrs.get("thetaMax")
        phiMin = dset.attrs.get("phiMin")
        phiMax = dset.attrs.get("phiMax")
        intensity = np.array( dset )

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    extent = np.array( [phiMin, phiMax, thetaMin, thetaMax] )
    extent = extent.reshape((-1,1))

    phi = np.linspace(phiMin, phiMax, intensity.shape[1])
    theta = np.linspace(thetaMin, thetaMax, intensity.shape[0])
    start = np.argmin( np.abs( phi-thetaMin) )
    end = np.argmin( np.abs( phi-thetaMax) )
    phi = phi[start:end]
    intensity = intensity[:,start:end]
    intensity = np.fliplr(intensity)
    im = ax.pcolor( phi, theta, intensity, cmap=cmap, norm=mpl.colors.LogNorm() )

    #ax.set_aspect( np.abs( (extent[1]-extent[0])/(extent[3]-extent[2]) ))
    ax.set_xlabel("Exit angle (deg)")
    ax.set_ylabel("Incident angle (deg)")
    fig.colorbar( im )
    fname = "Figures/incAngleSweep.jpeg"
    fig.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))

if __name__ == "__main__":
    main( sys.argv[1:] )
