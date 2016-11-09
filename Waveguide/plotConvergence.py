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
        error = np.array( hf.get("error") )
        nodes= np.array( hf.get("nodes") )

    slope, interscept, pvalue, rvalue, stderr = stats.linregress(np.log(nodes[1:]), np.log(error))
    print ("Convergence rate: %.2f"%(-slope))

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    #ax.plot( stepX/np.max(stepX), errorX, '^', color="black", label="X")
    #ax.plot( stepZ/np.max(stepZ), 0.1*errorZ, 's', color="black", label="Z")
    firstOrderConv = nodes.astype(np.float64)**(-1)
    firstOrderConv *= error[0]/firstOrderConv[1]
    secondOrderConv = nodes.astype(np.float64)**(-2)
    secondOrderConv *= error[0]/secondOrderConv[1]
    
    ax.plot( nodes[1:], error, 'o', fillstyle="none", color="black")
    ax.plot( nodes, firstOrderConv, color="black", ls="--" )
    ax.plot( nodes, secondOrderConv, color="black", ls="--")
    #ax.plot(stepRef, errorRef, color="black")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylabel("Error")
    ax.set_xlabel("Number of nodes")
    fname = "Figures/convergence.pdf"
    fig.savefig(fname, bbox_inches="tight")
    print ("Figure written to %s"%(fname))

if __name__ == "__main__":
    main(sys.argv[1:])
