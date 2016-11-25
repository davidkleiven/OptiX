import sys
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

def plot2D(data, stat, borders, field=None, phase=None):
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
        zmin = stat["zDiscretization"]["min"]*180.0/np.pi
        zmax = stat["zDiscretization"]["max"]*180.0/np.pi
        zlabel = "\$\\theta\$ (deg)"
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
    fig.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))
    if ( fname.find(".svg") != -1 ):
        psname = fname[:-3]+"pdf"
        subprocess.call(["inkscape", "--export-pdf=%s"%(psname), "--export-latex", fname])
        print ("PDF exported to %s"%(psname))


    plt.clf()
    frac = 1E-8
    maxval = np.max(np.abs(data)**2)
    maxval=1E-1
    minval = frac*maxval
    minval=1E-6
    #maxval = 1E-3 # T. Salditt et. al maxwav

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    im = ax.imshow(dataNorm, extent=extent, cmap=colormap, aspect=1.0, origin="lower", norm=mpl.colors.LogNorm(minval, maxval))
    ax.set_xlabel(zlabel)
    ax.set_ylabel(xlabel)
    ax.autoscale(False)
    fig.colorbar(im)
    if ( not borders is None ):
        borders.visualize( ax )
    ax.set_aspect( np.abs( (extent[1]-extent[0])/(extent[3]-extent[2]) ))
    #fname = "Figures/contourLogScale%d%s.jpeg"%(stat["UID"], appendName)
    fname = "Figures/contourLogScale%d%s.svg"%(stat["UID"], appendName)
    fig.savefig(fname, bbox_inches="tight")
    if ( fname.find(".svg") != -1 ):
        psname = fname[:-3]+"pdf"
        subprocess.call(["inkscape", "--export-pdf=%s"%(psname), "--export-latex", fname])
        print ("PDF exported to %s"%(psname))

    plt.show()
    print ("Figure written to %s"%(fname))

    plt.clf()
    if ( not field is None ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        im = ax.imshow(field, extent=extent, cmap=colormap, aspect=1.0, origin="lower")
        ax.set_xlabel(zlabel)
        ax.set_ylabel(xlabel)
        fig.colorbar( im )
        ax.autoscale(False)
        if ( not borders is None ):
            borders.visualize( ax )
        ax.set_aspect( np.abs( (extent[1]-extent[0])/(extent[3]-extent[2]) ))
        fname = "Figures/fieldLinScale%d.jpeg"%(stat["UID"])
        plt.show()
        fig.savefig(fname, bbox_inches="tight", dpi=800)
        print ("Figure written to %s"%(fname))

    if ( not phase is None ):
        #phase[np.abs(data)**2 < 1E-7] = np.NaN
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        im = ax.imshow(phase, extent=extent, cmap="viridis", aspect=1.0, origin="lower")
        ax.set_xlabel(zlabel)
        ax.set_ylabel(xlabel)
        fig.colorbar( im )
        ax.autoscale(False)
        if ( not borders is None ):
            borders.visualize( ax )
        ax.set_aspect( np.abs( (extent[1]-extent[0])/(extent[3]-extent[2]) ))
        fname = "Figures/phase%d.jpeg"%(stat["UID"])
        plt.savefig(fname, bbox_inches="tight", dpi=800)
        print ("Figure written to %s"%(fname))

def plot2Dsparse( x, z, intensity, stat ):
    print ("Using sparse plotting by triangulation...")
    #trianulation = tri.Triangulation( z/1000.0, x )
    zInterp = np.linspace(np.min(z), np.max(z), 501)
    xInterp = np.linspace(np.min(x), np.max(x), 501)
    Z, X = np.meshgrid( zInterp, xInterp )

    intensityInterp = interpolate.griddata( np.vstack((z,x)).T, intensity, (Z,X), method="linear")
    print ("Interpolation...")
    plt.clf()
    plt.contourf( Z/1000.0, X, intensityInterp**2, cmap="gist_heat")
    plt.xlabel("\$z\$ (\$\mathrm{\mu m}\$)")
    plt.ylabel("\$x\$ (nm)")
    plt.colorbar()
    fname = "Figures/contourLinScale.jpeg"
    fig.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))

    '''
    plt.clf()
    plt.tricontour( trianulation, intensity**2, cmap="gist_heat", norm=mpl.colors.LogNorm(np.min(intensity**2), np.max(intensity**2)))
    plt.xlabel("\$z\$ (\$\mathrm{\mu m}\$)")
    plt.ylabel("\$x\$ (nm)")
    plt.colorbar()
    fname = "Figures/contourLogScale.jpeg"
    plt.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))
    '''

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
            data = np.array( hf.get("dataset") )

    fieldData = None
    try:
        with h5.File(stat["fieldData"], 'r') as hf:
            fieldData = np.array( hf.get("dataset") )
    except:
        fieldData = None

    phaseData = None
    try:
        with h5.File(stat["phase"], 'r') as hf:
            phaseData = np.array( hf.get("dataset") )
    except Exception as exc:
        print str(exc)
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

    '''
    if ( np.min(xInside) > stat["xDiscretization"]["min"]):
        x0 = stat["xDiscretization"]["min"]
    else:
        x0 = np.min(xInside)
    stat["x0"] = x0
    '''
    stat["x0"] = stat["xDiscretization"]["min"]

    if ( stat["sparseSave"] ):
        plot2Dsparse( xVal, zVal, intensity, stat )
    else:
        data = data.T # Transpose the dataset
        fieldData = fieldData.T
        if ( not phaseData is None ):
            phaseData = phaseData.T
        plot2D( data, stat, borders, field=fieldData, phase=phaseData )
    #plotWG( xInside-x0, zInside )

    # Plot transmission. Put in try catch as some of the simulaitons do not compute the transmission
    try:
        with h5.File(stat["Transmission"]["file"], 'r') as hf:
            data = np.array( hf.get("transmission") )
        trans.plotTransmission( data, stat )
    except Exception as exc:
        print str(exc)
        print ("Error when plotting transmission")

if __name__ == "__main__":
    main(sys.argv[1:])
