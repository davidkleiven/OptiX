import sys
sys.path.append("../FresnelFDTD")
sys.path.append("../")
import mplLaTeX as ml
import matplotlib as mpl
#mpl.rcParams.update(ml.params)
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 28
import numpy as np
import json
import h5py as h5
from matplotlib import pyplot as plt
import subprocess
import colorScheme as cs
try:
    import colormaps as cmaps
    cmap = cmaps.viridis
except:
    cmap = "viridis"

#cmap="jet"
DELTA = 4.14E-5 # Salditt et al
BETA = 3.45E-6 # Salditt et al

def main( argv ):
    global cmap
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        elif ( arg.find("--help") != -1 ):
            print ("Usage: python incAngleSweep.py --file=<hdf5 datafile>")
            print ("file: HDF5 file created by the incidentAngleSweep application")
            print ("help: Print this message")
            print ("cmap: Specify the cmap (default: viridis)")
            return
        elif ( arg.find("--cmap=") != -1 ):
            cmap = arg.split("--cmap=")[1]
        else:
            print ("Unknown argument %s"%(arg))
            return

    with h5.File( fname, 'r' ) as hf:
        dset = hf.get("intensity")
        thetaMin = float( dset.attrs.get("thetaMin") )
        thetaMax = float( dset.attrs.get("thetaMax") )
        phiMin = float( dset.attrs.get("phiMin") )
        phiMax = float( dset.attrs.get("phiMax") )
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
    #intensity = np.fliplr(intensity)
    extent = [phiMin, phiMax, thetaMin, thetaMax]
    intensity = intensity.T
    im = ax.imshow( intensity, extent=extent, cmap=cmap, norm=mpl.colors.LogNorm(), origin="lower" )
    ax.set_aspect( (extent[1]-extent[0])/(extent[3]-extent[2]))

    #ax.set_aspect( np.abs( (extent[1]-extent[0])/(extent[3]-extent[2]) ))
    ax.set_xlabel("Exit angle (deg)")
    ax.set_ylabel("Incident angle (deg)")
    ax.autoscale(False)
    fig.colorbar( im )
    fname = "Figures/incAngleSweep.svg"
    fig.savefig(fname, bbox_inches="tight", dpi=800)
    psname = "Figures/incAngleSweep.ps"
    subprocess.call(["inkscape", "--export-ps=%s"%(psname), "--export-latex", fname])
    print ("Figure written to %s"%(fname))
    plt.show()

if __name__ == "__main__":
    main( sys.argv[1:] )
