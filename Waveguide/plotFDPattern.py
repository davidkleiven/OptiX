import sys
sys.path.append("../FresnelFDTD")
import mplLaTeX as ml
import matplotlib as mpl
mpl.rcParams.update(ml.params)
import numpy as np
import h5py as h5
import json
from matplotlib import pyplot as plt

def plot2D(data, stat):
    x = np.linspace(stat["xDiscretization"]["min"], stat["xDiscretization"]["max"], data.shape[0])
    z = np.linspace(stat["zDiscretization"]["min"], stat["zDiscretization"]["max"], data.shape[1])
    Z,X = np.meshgrid(z,x)

    print np.max(data)
    print np.min(data)

    for i in range(0, data.shape[0]):
        print data[i,:]

    #plt.contourf(Z/1000.0,X,data, 200, cmap="gist_heat", norm=mpl.colors.LogNorm())
    plt.contourf(data, cmap="gist_heat", norm=mpl.colors.LogNorm())
    plt.xlabel("$z$ ($\mathrm{\mu m}$)")
    plt.ylabel("$x$ (nm)")
    fname = "Figures/contourLinScale.jpeg"
    plt.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))

def main(argv):
    fname = ""
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        else:
            print ("Unknown argument %s"%(arg))
            return 1

    if ( fname == ""):
        print ("No json file specified")
        return 1

    try:
        infile = open(fname, "r")
        stat = json.load(infile)
    except:
        print("Could not open file %s"%(fname))
        return 1

    with h5.File(stat["datafile"], "r") as hf:
        dataReal = np.array( hf.get("real") )
        dataImag = np.array( hf.get("imag") )
    print dataReal.shape
    intensity = dataReal**2 + dataImag**2
    plot2D( intensity, stat )

if __name__ == "__main__":
    main(sys.argv[1:])
