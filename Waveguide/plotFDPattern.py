import sys
from PLOD import controlGUI as cg
import tkinter as tk
import subprocess
sys.path.append("../FresnelFDTD")
import mplLaTeX as ml
import matplotlib as mpl
#mpl.rcParams.update(ml.params)
mpl.rcParams["svg.fonttype"]="none"
mpl.rcParams["axes.linewidth"] = 0.1
mpl.rcParams["font.size"] = 28
CALL_OVER_SSH = False
if ( CALL_OVER_SSH ):
    mpl.use("Agg")
import numpy as np
import h5py as h5
import json
from matplotlib import pyplot as plt
from scipy import interpolate
import transmission as trans
import waveguideBorder as wgb

def plot2D(data, xmin, xmax, zmin, zmax, uid, borders, control, field=None, phase=None):
    colormap="viridis"
    colormap="nipy_spectral"

    zmin /= 1E6
    zmax /= 1E6

    zlabel = "\$ v (\SI{}{\milli\meter}\)$"
    xlabel = "\$u (\SI{}{\\nano\meter})\$"

    normalize = True
    appendName = ""
    if ( normalize ):
        appendName = "normalized"

    extent = [zmin, zmax, xmin, xmax]
    dataNorm = np.abs(data)**2
    if ( normalize ):
        dataNorm /= np.max(dataNorm)

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    im = ax.imshow(dataNorm, extent=extent, cmap=colormap, aspect=1.0, origin ="lower")
    ax.set_xlabel(zlabel)
    ax.set_ylabel(xlabel)
    ax.autoscale(False)
    ax.set_aspect( np.abs( (extent[1]-extent[0])/(extent[3]-extent[2]) ))
    fig.colorbar( im )
    if ( not borders is None ):
        borders.visualize( ax )
    #fname = "Figures/contourLinScale%d%s.jpeg"%(stat["UID"], appendName)
    fname = "Figures/contourLinScale%d%s.svg"%(uid, appendName)
    control.attach( fig, ax, fname )

    frac = 1E-8
    maxval = np.max(np.abs(data)**2)
    #maxval=1E-1
    minval = frac*maxval
    #minval=1E-6
    #maxval = 1E-3 # T. Salditt et. al maxwav

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(1,1,1)
    im = ax2.imshow(dataNorm, extent=extent, cmap=colormap, aspect=1.0, origin="lower", norm=mpl.colors.LogNorm(minval, maxval))
    ax2.set_xlabel(zlabel)
    ax2.set_ylabel(xlabel)
    ax2.autoscale(False)
    fig2.colorbar(im)
    if ( not borders is None ):
        borders.visualize( ax )
    ax2.set_aspect( np.abs( (extent[1]-extent[0])/(extent[3]-extent[2]) ))
    fname = "Figures/contourLogScale%d%s.svg"%(uid, appendName)
    control.attach( fig2, ax2, fname )

    if ( not field is None ):
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(1,1,1)
        im = ax3.imshow(field, extent=extent, cmap=colormap, aspect=1.0, origin="lower")
        ax3.set_xlabel(zlabel)
        ax3.set_ylabel(xlabel)
        fig3.colorbar( im )
        ax3.autoscale(False)
        if ( not borders is None ):
            borders.visualize( ax )
        ax3.set_aspect( np.abs( (extent[1]-extent[0])/(extent[3]-extent[2]) ))
        fname = "Figures/fieldLinScale%d.jpeg"%(uid)
        control.attach( fig3, ax3, fname )

    if ( not phase is None ):
        #phase[np.abs(data)**2 < 1E-7] = np.NaN
        fig4 = plt.figure()
        ax4 = fig4.add_subplot(1,1,1)
        im = ax4.imshow(phase, extent=extent, cmap=colormap, aspect=1.0, origin="lower")
        ax4.set_xlabel(zlabel)
        ax4.set_ylabel(xlabel)
        fig4.colorbar( im )
        ax4.autoscale(False)
        if ( not borders is None ):
            borders.visualize( ax )
        ax4.set_aspect( np.abs( (extent[1]-extent[0])/(extent[3]-extent[2]) ))
        fname = "Figures/phase%d.jpeg"%(uid)
        control.attach( fig4, ax4, fname )

def plotWG( x, z ):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot( z/1000.0, x,'k.')
    ax.set_xlabel("$z$ $\mathrm{\mu m}$")
    ax.set_ylabel("$x$ (nm)")
    fname = "Figures/pointsInWG.jpeg"
    fig.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))

def readBorders( hf, crdsyst ):
    borders = wgb.WaveGuideBorders()
    mxBorders = 10
    if ( crdsyst == "cylindrical" ):
        factor = 180.0/np.pi
    else:
        factor = 1E-3
    for i in range(0, mxBorders):
        uxname = "upperBorderX%d"%(i)
        uzname = "upperBorderZ%d"%(i)
        bxname = "lowerBorderX%d"%(i)
        bzname = "lowerBorderZ%d"%(i)
        x1 = hf.get(bxname)
        z1 = hf.get(bzname)
        x2 = hf.get(uxname)
        z2 = hf.get(uzname)
        if ( x1 is None ) or ( z1 is None ) or ( x2 is None ) or ( z2 is None ):
            print ("Read %d waveguide borders"%(i))
            return borders

        borders.addBorder(np.array(x1), np.array(z1)*factor, np.array(x2), np.array(z2)*factor)
    return borders

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

    root = tk.Tk()
    control = cg.Control( root )

    zmin = 0.0
    zmax = 0.0
    uid = 0
    with h5.File( fname, "r") as hf:
        group = hf.get("/data")
        data = np.array( hf.get("/data/amplitude") )
        if ( "/data/phase" in hf.keys() ):
            phaseData = np.array( hf.get("/data/phase") )
        else:
            phaseData = None
        if ( "/data/transmittivity" in hf.keys() ):
            transData = np.array( hf.get("/data/transmittivity") )
        else:
            transData = None
        zmin = group.attrs.get("zmin")
        zmax = group.attrs.get("zmax")
        xmin = group.attrs.get("xmin")
        xmax = group.attrs.get("xmax")
        uid = group.attrs.get("uid")

    # TODO: Fix the borders!
    borders = None
    try:
        crdsyst = stat["crd"]
    except:
        crdsyst = "cartesian"
    try:
        with h5.File(stat["wgfile"], 'r') as hf:
            borders = readBorders( hf, crdsyst )
    except:
        print ("Borders were not found!")

    plot2D( data, xmin, xmax, zmin, zmax, uid, borders, control, field=None, phase=phaseData )
    if ( not transData is None ):
        # Due to downsampling and merging of some datasets there can be some elements that are not set
        # Quick fix: Just remove, say, the last 5 elements of the array
        transData = transData[:-5]
        trans.plotTransmission( transData, zmin, zmax, uid, control )

    root.mainloop()

if __name__ == "__main__":
    main(sys.argv[1:])
