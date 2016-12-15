from __future__ import print_function
import numpy as np
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 28
mpl.rcParams["axes.unicode_minus"] = False
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
    setLabels = True
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3"]
    for entry in data["dsets"]:
        simulations = {}
        for pos in positions:
            simulations[str(int(pos))] = Data()

        for uid in entry["uids"]:
            fname = data["basename"]+str(uid)+".json"
            infile = open( fname, 'r')
            ctl = json.load( infile )
            infile.close()

            with h5.File( ctl["Transmission"]["file"], 'r' ) as hf:
                T = np.array( hf.get("transmission") )
            z = np.linspace( ctl["zDiscretization"]["min"], ctl["zDiscretization"]["max"], len(T) )/1000.0
            if ( ctl["waveguide"]["crd"] == "cylindrical" ):
                z *= ctl["waveguide"]["RadiusOfCurvature"]

            for pos in positions:
                indx = np.argmin( np.abs(z-float(pos)) )
                simulations[str(int(pos))].dx.append(ctl["xDiscretization"]["step"])
                simulations[str(int(pos))].T.append(T[indx])
        indx = 0
        for pos in positions:
            dx = np.array( simulations[str(int(pos))].dx )
            T = np.array( simulations[str(int(pos))].T )
            color = colors[indx%len(colors)]
            if ( setLabels ):
                ax.plot( 1.0/dx, np.log(T), label="%.2f mm"%(pos/1000.0), color=color, ls=entry["ls"], marker=entry["marker"], ms=7, fillstyle="none" )
            else:
                ax.plot( 1.0/dx, np.log(T), ls=entry["ls"], color=color, marker=entry["marker"], ms=7, fillstyle="none" )
            indx += 1
        setLabels = False
    ax.legend(bbox_to_anchor=[1.05,0.45], frameon=False, labelspacing=0.005)
    ax.set_xlabel("\$\Delta x^{-1}\$ (nm)\$^{-1}\$")
    ax.set_ylabel("\$\log T\$")
    plt.show()

    # Save figure
    figname = data["figname"]
    fig.savefig(figname)
    print ("Figure written to %s"%(figname))
    if ( figname[-3:] == "svg" ):
        psname = figname[:-3]+"ps"
        subprocess.call( ["inkscape", "--export-ps=%s"%(psname), "--export-latex", figname] )
        print("Figure written to %s"%(psname))
if __name__ == "__main__":
    main( sys.argv[1:] )
