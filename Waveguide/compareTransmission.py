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
    MSG = "Usage: python compareTransmission.py --file=<jsonfile> [--help]\n"
    MSG += "help: Print this message\n"
    MSG += "file: Json files describing the UIDs to plot"
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        elif ( arg.find("--help") != -1 ):
            print MSG
            return 0
        else:
            print ("Unknown argument %s"%(arg))
            return 0

    try:
        infile = open(fname,'r')
        param = json.load(infile)
        infile.close()
    except Exception as exc:
        print str(exc)
        print ("Error when opening/parsing file %s"%(fname))
        return 0

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    indx = 0
    markers = ["s", "^", '.']
    minOfMaxZ = np.inf
    for entry in param["entries"]:
        ctlfile = param["basename"]+"%d.json"%(entry["uid"])
        try:
            infile = open(ctlfile, 'r')
            stat = json.load(infile)
            infile.close()
        except Exception as exc:
            print str(exc)
            print ("Error when opening/parsing file %s"%(ctlfile))
            return 0

        with h5.File(stat["Transmission"]["file"], 'r') as hf:
            data = np.array( hf.get( "transmission" ) )

        #data = data[::param["numberOfPoints"]]
        z = np.linspace(stat["Transmission"]["zStart"], stat["Transmission"]["zEnd"], len(data))
        if ( stat["Transmission"]["zEnd"] < minOfMaxZ ):
            minOfMaxZ = stat["Transmission"]["zEnd"]/1E6
        ax.plot( z/1E6, np.log(data), color="black", marker=markers[indx], ms=3, linestyle="None", label=entry["label"])
        indx += 1

    ax.set_xlabel("$z$ (mm)")
    ax.set_xlim(right=minOfMaxZ)
    ax.set_ylabel("$\ln$ (Transmission)")
    ax.legend(loc="upper right", frameon=False)
    fname = "Figures/"+param["figurename"]
    fig.savefig(fname, bbox_inches="tight")
    print ("Figure written to %s"%(fname))

if __name__ == "__main__":
    main(sys.argv[1:])
