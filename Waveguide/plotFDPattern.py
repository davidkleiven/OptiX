import sys
sys.path.append("../FresnelFDTD")
import mplLaTeX as ml
import matplotlib as mpl
mpl.rcParams.update(ml.params)
import numpy as np
import h5py as h5
import json
from matplotlib import pyplot as plt

def plot2D(dataReal, dataImag, stat):
    x = np.linspace(stat["xDiscretization"]["min"], stat["xDiscretization"]["max"], dataReal.shape[0])
    x -= stat["x0"]
    z = np.linspace(stat["zDiscretization"]["min"], stat["zDiscretization"]["max"], dataReal.shape[1])
    Z,X = np.meshgrid(z,x)

    k = 2.0*np.pi/0.1569
    field = (dataReal+1j*dataImag)
    plt.contourf(Z/1000.0,X, np.abs(field)**2, 200, cmap="gist_heat")
    plt.xlabel("$z$ ($\mathrm{\mu m}$)")
    plt.ylabel("$x$ (nm)")
    plt.colorbar()
    fname = "Figures/contourLinScale.jpeg"
    plt.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))

    plt.clf()
    plt.contourf(Z/1000.0,X, np.abs(field)**2, 200, cmap="gist_heat", norm=mpl.colors.LogNorm())
    plt.xlabel("$z$ ($\mathrm{\mu m}$)")
    plt.ylabel("$x$ (nm)")
    plt.colorbar()
    fname = "Figures/contourLogScale.jpeg"
    plt.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))

def plotWG( x, z ):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot( z/1000.0, x, 'k.')
    ax.set_xlabel("$z$ $\mathrm{\mu m}$")
    ax.set_ylabel("$x$ (nm)")
    fname = "Figures/pointsInWG.jpeg"
    fig.savefig(fname, bbox_inches="tight", dpi=800)
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

    with h5.File(stat["wgfile"], 'r') as hf:
        xInside = np.array( hf.get("xInside"))
        zInside = np.array( hf.get("zInside"))

    if ( np.min(xInside) > stat["xDiscretization"]["min"]):
        x0 = stat["xDiscretization"]["min"]
    else:
        x0 = np.min(xInside)
    stat["x0"] = x0
    plot2D( dataReal, dataImag, stat )
    plotWG( xInside-x0, zInside )

if __name__ == "__main__":
    main(sys.argv[1:])
