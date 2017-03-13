import numpy as np
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 18
mpl.rcParams["axes.unicode_minus"]=False
from matplotlib import pyplot as plt
import h5py as h5

def cleanFileWebPlotDig( x, y ):
    xnew = np.zeros(len(x))
    ynew = np.zeros(len(y))
    nInserted = 0
    for i in range(0,len(x)-1):
        if ( x[i+1] > x[i] ):
            xnew[nInserted] = x[i]
            ynew[nInserted] = y[i]
            nInserted += 1
    return xnew[:nInserted], ynew[:nInserted]

def main():
    fname = "data/sphereCrossSection.h5"
    fnameWiki = "data/wikipediaRCSMetallicSphere.csv"
    data = np.loadtxt( fnameWiki, delimiter="," )
    with h5.File(fname, 'r') as hf:
        refPlane = np.array( hf.get("refPlaneFlux") )
        boxFluxScat = np.array( hf.get("scatFluxBox") )
        freq = np.array( hf.get("freq") )
        eps = np.array( hf.get("eps") )

    R = 10.0
    f = np.linspace( freq[0], freq[1]*freq[2], freq[1])
    k = 2.0*np.pi*f
    ratio = np.abs( boxFluxScat/refPlane )

    exactX, exactY = cleanFileWebPlotDig( data[:,0], data[:,1] )
    # Normalize
    ratio *= np.sum(data[:,1])*len(ratio)/(np.sum(ratio)*len(data[:,1]))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot( k*R, ratio, color="#e41a1c", label="Simulated")
    ax.plot( exactX, exactY, color="#377eb8", label="Mie")
    ax.legend(loc="lower left", frameon=False, labelspacing=0.05)
    ax.set_yscale("log")
    ax.set_xlabel("kR")
    ax.set_ylabel("Scattering cross section")
    ax.set_xscale("log")

    # Plot projections
    fig2 = plt.figure()
    proj = np.sum(eps, axis=2)
    ax1 = fig2.add_subplot(2,2,1)
    ax1.imshow(proj, cmap="bone")
    proj = np.sum(eps,axis=1)
    ax2 = fig2.add_subplot(2,2,2)
    ax2.imshow(proj, cmap="bone")
    ax3 = fig2.add_subplot(2,2,3)
    proj = np.sum(eps,axis=0)
    ax3.imshow(proj, cmap="bone")
    plt.show()

if __name__ == "__main__":
    main()
