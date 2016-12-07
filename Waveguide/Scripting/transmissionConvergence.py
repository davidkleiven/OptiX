from __future__ import print_function
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import h5py as h5
import subprocess
import sys
import json

try:
    import colormaps as cmaps
    cmap = cmaps.viridis
except:
    cmap = "viridis"

class Data:
    def __init__(self):
        self.dx = []
        self.T = []

def main( argv ):
    fname = ""
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]

    if ( fname == "" ):
        print ("No filename given!")
        return

    infile = open( fname, 'r' )
    data = json.load( infile )
    infile.close()

    positions = data["positions"]
    values = []
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for entry in data["dsets"]:
        simulations = {}
        for pos in positions:
            simulations[str(int(pos))] = Data()

        for uid in entry["uids"]:
            fname = data["basename"]+entry["uid"]+".json"
            infile = open( fname, 'r')
            ctl = json.load( infile )
            infile.close()

            with h5.File( ctl["Transmission"]["file"], 'r' ) as hf:
                T = np.array( hf.get("transmission") )
            z = np.linspace( ctl["zDiscretization"]["min"], ctl["zDiscretization"]["max"], len(T) )
            for pos in positions:
                indx = np.argmin( np.abs(z-float(pos)) )
                simulations[str(int(pos))].dx.append(ctl["xDiscretization"]["step"])
                simulations[str(int(pos))].T.append(T[indx])
        for pos in positions:
            dx = np.array( simulations[str(int(pos))].dx )
            T = np.array( simulations[str(int(pos))].T )
            ax.plot( 1.0/dx, T )
    plt.show()

if __name__ == "__main__":
    main( sys.argv[1:] )
