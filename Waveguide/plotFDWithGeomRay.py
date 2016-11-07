import sys
sys.path.append("../FresnelFDTD")
sys.path.append("GeometricalOptics")
import mplLaTeX as ml
import matplotlib as mpl
#mpl.rcParams.update(ml.params)
import numpy as np
import h5py as h5
import json
from matplotlib import pyplot as plt
from scipy import interpolate
import transmission as trans
import waveguideBorder as wgb
import circularWGCartesian as cwgc
import plotFDPattern as pfdp

def main(argv):
    fname = ""
    colormap = "viridis"
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

    fieldData = None
    try:
        with h5.File(stat["fieldData"], 'r') as hf:
            fieldData = np.array( hf.get("dataset") )
    except:
        fieldData = None
        return 1

    # Start plotting
    try:
        geom = cwgc.CurvedWGCartesian()
        geom.R = stat["waveguide"]["RadiusOfCurvature"]
        geom.width = stat["waveguide"]["Width"]
    except:
        print ("Error when reading wavguide parameters")
        return 1

    # Read borders
    with h5.File(stat["wgfile"], 'r') as hf:
        borders = pfdp.readBorders( hf, "cartesian" )

    extent = [stat["zDiscretization"]["min"], stat["zDiscretization"]["max"], stat["xDiscretization"]["min"], stat["xDiscretization"]["max"]]
    k = 2.0*np.pi/0.1569
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    im = ax.imshow(np.abs(fieldData)**2, extent=extent, cmap=colormap, origin ="lower")
    ax.set_xlabel(zlabel)
    ax.set_ylabel(xlabel)
    ax.set_aspect( np.abs( (extent[1]-extent[0])/(extent[3]-extent[2]) ))
    fig.colorbar( im )

    borders.visualize( ax )
    geom.solve( 50.0, 0.0, [0.0,1.0])
    geom.plot( ax )
    fname = "Figures/contourLinScaleWithGeomRay%d.jpeg"%(stat["UID"])
    fig.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))

if __name__ == "__main__":
    main(sys.argv[1:])
