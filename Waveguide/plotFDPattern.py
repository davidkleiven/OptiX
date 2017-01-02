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

def plot2D(data, stat, borders, control, field=None, phase=None):
    colormap="viridis"
    print ("Plotting the full matrix...")
    '''
    x = np.linspace(stat["xDiscretization"]["min"], stat["xDiscretization"]["max"], data.shape[0])
    x -= stat["x0"]
    z = np.linspace(stat["zDiscretization"]["min"], stat["zDiscretization"]["max"], data.shape[1])
    Z,X = np.meshgrid(z,x)
    '''

    normalize = True
    appendName = ""
    if ( normalize ):
        appendName = "normalized"

    try:
        crd = stat["waveguide"]["crd"]
    except:
        crd="cartesian"

    if ( crd == "cylindrical" ):
        zmin = stat["zDiscretization"]["min"]*stat["waveguide"]["RadiusOfCurvature"]/1000.0
        zmax = stat["zDiscretization"]["max"]*stat["waveguide"]["RadiusOfCurvature"]/1000.0
        zlabel = "\$R\\theta (\SI{}{\mircro\meter})\$"
        xlabel = "\$r\$ (nm)"
    else:
        zmin = stat["zDiscretization"]["min"]/1000.0
        zmax = stat["zDiscretization"]["max"]/1000.0
        zlabel = "\$z\$ (\$\micro\$m)"
        xlabel = "\$x\$ (nm)"

    extent = [zmin, zmax, stat["xDiscretization"]["min"], stat["xDiscretization"]["max"]]
    dataNorm = np.abs(data)**2
    if ( normalize ):
        dataNorm /= np.max(dataNorm)

    k = 2.0*np.pi/0.1569
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
    fname = "Figures/contourLinScale%d%s.svg"%(stat["UID"], appendName)
    control.attach( fig, ax, fname )

    frac = 1E-8
    maxval = np.max(np.abs(data)**2)
    maxval=1E-1
    minval = frac*maxval
    minval=1E-6
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
    #fname = "Figures/contourLogScale%d%s.jpeg"%(stat["UID"], appendName)
    fname = "Figures/contourLogScale%d%s.svg"%(stat["UID"], appendName)
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
        fname = "Figures/fieldLinScale%d.jpeg"%(stat["UID"])
        control.attach( fig3, ax3, fname )

    if ( not phase is None ):
        #phase[np.abs(data)**2 < 1E-7] = np.NaN
        fig4 = plt.figure()
        ax4 = fig4.add_subplot(1,1,1)
        im = ax4.imshow(phase, extent=extent, cmap="viridis", aspect=1.0, origin="lower")
        ax4.set_xlabel(zlabel)
        ax4.set_ylabel(xlabel)
        fig4.colorbar( im )
        ax4.autoscale(False)
        if ( not borders is None ):
            borders.visualize( ax )
        ax4.set_aspect( np.abs( (extent[1]-extent[0])/(extent[3]-extent[2]) ))
        fname = "Figures/phase%d.jpeg"%(stat["UID"])
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
    try:
        infile = open(fname, "r")
        stat = json.load(infile)
    except:
        print("Could not open file %s"%(fname))
        return 1

    try:
        val = stat["sparseSave"]
    except:
        stat["sparseSave"] = False

    if ( stat["sparseSave"] ):
        with h5.File(stat["datafile"], 'r') as hf:
            xVal = np.array( hf.get("x") )
            zVal = np.array( hf.get("z") )
            intensity = np.array( hf.get("intensity") )
    else:
        with h5.File(stat["datafile"], "r") as hf:
            data = np.array( hf.get("amplitude") )

    fieldData = None
    try:
        with h5.File(stat["datafile"], 'r') as hf:
            fieldData = np.array( hf.get("amplitude") )
    except:
        fieldData = None

    phaseData = None
    try:
        with h5.File(stat["datafile"], 'r') as hf:
            phaseData = np.array( hf.get("phase") )
    except Exception as exc:
        print ( str(exc) )
        phaseData = None

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


    stat["x0"] = stat["xDiscretization"]["min"]

    plot2D( data, stat, borders, control, field=fieldData, phase=phaseData )

    # Plot transmission. Put in try catch as some of the simulaitons do not compute the transmission
    try:
        with h5.File(stat["datafile"], 'r') as hf:
            dset = hf.get("transmittivity")
            data = np.array( dset )
            zmin = float( dset.attrs.get("zmin") )
            zmax = float( dset.attrs.get("zmax") )
            uid = int( dset.attrs.get("uid") )
        trans.plotTransmission( data, zmin, zmax, uid, control )
    except Exception as exc:
        print ( str(exc) )
        print ("Error when plotting transmission")
    root.mainloop()

if __name__ == "__main__":
    main(sys.argv[1:])
