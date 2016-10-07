import sys
sys.path.append("../FresnelFDTD")
import mplLaTeX as ml
import matplotlib as mpl
mpl.rcParams.update(ml.params)
import numpy as np
import json
import h5py as h5
from matplotlib import pyplot as plt

def main(argv):
    fname = ""
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        else:
            print ("Unknown argument %s"%(fname))
            return

    if ( fname == "" ):
        print ("No file specified")
        return

    infile = open( fname , "r" )
    stat = json.load(infile)
    infile.close()

    with h5.File(stat["solutionfile"], "r") as hf:
        data = np.array( hf.get(stat["datasetname"]) )

    u = np.linspace( stat["solver"]["xmin"], stat["solver"]["xmax"], len(data))

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot( u, data, 'k')
    figname = "Figures/profile.pdf"
    fig.savefig(figname, bbox_inches="tight")
    print ("Figure written to %s"%(figname))

    # Plot potential
    with h5.File(stat["potentialFname"], "r") as hf:
        potential = np.array( hf.get(stat["potentialLabel"]))

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(potential, 'k')
    fname = "Figures/potential.pdf"
    fig.savefig( fname, bbox_inches="tight")
    print ("Figure written to %s"%(fname))
if __name__ == "__main__":
    main(sys.argv[1:])
