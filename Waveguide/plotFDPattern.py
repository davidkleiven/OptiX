import sys
sys.path.append("../FresnelFDTD")
import mplLaTeX as ml
import matplotlib as mpl
#mpl.rcParams.update(ml.params)
CALL_OVER_SSH = True
if ( CALL_OVER_SSH ):
    mpl.use("Agg")
import numpy as np
import h5py as h5
import json
from matplotlib import pyplot as plt
from scipy import interpolate
import transmission as trans

def plot2D(data, stat, field=None):
    print ("Plotting the full matrix...")
    x = np.linspace(stat["xDiscretization"]["min"], stat["xDiscretization"]["max"], data.shape[0])
    x -= stat["x0"]
    z = np.linspace(stat["zDiscretization"]["min"], stat["zDiscretization"]["max"], data.shape[1])
    Z,X = np.meshgrid(z,x)
    extent = [stat["zDiscretization"]["min"]/1000.0, stat["zDiscretization"]["max"]/1000.0,
                 stat["xDiscretization"]["min"], stat["xDiscretization"]["max"]]

    k = 2.0*np.pi/0.1569
    plt.imshow(np.abs(data)**2, extent=extent, cmap="coolwarm", aspect=1.0, origin ="lower")
    plt.xlabel("$z$ ($\mathrm{\mu m}$)")
    plt.ylabel("$x$ (nm)")
    plt.gca().set_aspect( np.abs( (extent[1]-extent[0])/(extent[3]-extent[2]) ))
    plt.colorbar()
    fname = "Figures/contourLinScale%d.jpeg"%(stat["UID"])
    plt.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))

    plt.clf()
    minval = 1E-8
    maxval = np.max(np.abs(data)**2)
    maxval = 1E-3 # T. Salditt et. al maxwav

    plt.imshow(np.abs(data)**2, extent=extent, cmap="coolwarm", aspect=1.0, origin="lower", norm=mpl.colors.LogNorm(minval, maxval))
    plt.xlabel("$z$ ($\mathrm{\mu m}$)")
    plt.ylabel("$x$ (nm)")
    plt.colorbar()
    plt.gca().set_aspect( np.abs( (extent[1]-extent[0])/(extent[3]-extent[2]) ))
    fname = "Figures/contourLogScale%d.jpeg"%(stat["UID"])
    plt.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))

    plt.clf()
    if ( not field is None ):
        plt.imshow(field, extent=extent, cmap="coolwarm", aspect=1.0, origin="lower")
        plt.xlabel("$z$ ($\mathrm{\mu m}$)")
        plt.ylabel("$x$ (nm)")
        plt.colorbar()
        plt.gca().set_aspect( np.abs( (extent[1]-extent[0])/(extent[3]-extent[2]) ))
        fname = "Figures/fieldLinScale%d.jpeg"%(stat["UID"])
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
    plt.xlabel("$z$ ($\mathrm{\mu m}$)")
    plt.ylabel("$x$ (nm)")
    plt.colorbar()
    fname = "Figures/contourLinScale.jpeg"
    plt.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))

    '''
    plt.clf()
    plt.tricontour( trianulation, intensity**2, cmap="gist_heat", norm=mpl.colors.LogNorm(np.min(intensity**2), np.max(intensity**2)))
    plt.xlabel("$z$ ($\mathrm{\mu m}$)")
    plt.ylabel("$x$ (nm)")
    plt.colorbar()
    fname = "Figures/contourLogScale.jpeg"
    plt.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))
    '''

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

    with h5.File(stat["wgfile"], 'r') as hf:
        xInside = np.array( hf.get("xInside"))
        zInside = np.array( hf.get("zInside"))

    if ( np.min(xInside) > stat["xDiscretization"]["min"]):
        x0 = stat["xDiscretization"]["min"]
    else:
        x0 = np.min(xInside)
    stat["x0"] = x0

    if ( stat["sparseSave"] ):
        plot2Dsparse( xVal, zVal, intensity, stat )
    else:
        data = data.T # Transpose the dataset
        fieldData = fieldData.T
        plot2D( data, stat, field=fieldData )
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
