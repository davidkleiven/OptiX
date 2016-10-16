import sys
sys.path.append("../FresnelFDTD")
import mplLaTeX as ml
import matplotlib as mpl
mpl.rcParams.update(ml.params)
import numpy as np
import h5py as h5
import json
from matplotlib import pyplot as plt
from scipy import stats

def main( argv ):
    MSG = "Usage: python plotConvergence.py --file=<hdf5 file> [--help]\n"
    MSG += "help: Print this message\n"
    MSG += "file: HDF5 file containingn the convergence data"
    for arg in argv:
        if ( arg.find("--help") != -1 ):
            print MSG
            return 0
        elif ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        else:
            print ("Unknown argument %s"%(arg))
            return 1

    with h5.File(fname, 'r') as hf:
        stepX = np.array( hf.get("stepsizeX") )
        errorX = np.array( hf.get("errorRatioX") )
        stepZ = np.array( hf.get("stepsizeZ") )
        errorZ = np.array( hf.get("errorRatioZ") )

    stepsize = np.zeros(len(errorX))
    for i in range(0, len(stepsize)):
        stepsize[i] = 2**i
    stepRef = np.linspace(0.01, 1.0, 101)
    errorRef = 10*stepRef**2
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    #ax.plot( stepX/np.max(stepX), errorX, '^', color="black", label="X")
    #ax.plot( stepZ/np.max(stepZ), 0.1*errorZ, 's', color="black", label="Z")
    ax.plot( stepsize, errorX, '^', color="black", label="X")
    ax.plot( stepsize, errorZ, 's', color="black", label="Z")
    #ax.plot(stepRef, errorRef, color="black")
    ax.set_yscale("log")
    ax.set_xscale("log", basex=2)
    ax.set_ylabel("Error")
    ax.set_xlabel("Normalized stepsize")
    ax.legend(loc="lower right", frameon=False)
    fname = "Figures/convergence.pdf"
    fig.savefig(fname, bbox_inches="tight")
    print ("Figure written to %s"%(fname))

if __name__ == "__main__":
    main(sys.argv[1:])
