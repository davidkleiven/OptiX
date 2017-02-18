import sys
import tkinter as tk
import proj3D
import h5py as h5
from PLOD import controlGUI as cg
import numpy as np

def main( argv ):
    fname = ""
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        else:
            print ("Unrecognized option %s"%(arg))
            return

    if ( fname == "" ):
        print ("No file given!")
        return

    lim = proj3D.Limits()
    ampPlot = proj3D.Plotter3D()
    phasePlot = proj3D.Phase()
    ffPlot = proj3D.FarField()
    intensityPlot = proj3D.IntensityPlot()

    with h5.File( fname, 'r' ) as hf:
        group = hf.get("/data")
        ffset = hf.get("/data/farField")
        ffPlot.data = np.array( ffset ).T
        phasePlot.data = np.array( hf.get("/data/exitPhase") ).T
        intensityPlot.data = np.array( hf.get("/data/exitIntensity") ).T
        lim.zmin = group.attrs.get("zmin")
        lim.zmax = group.attrs.get("zmax")
        lim.xmin = group.attrs.get("xmin")/1000.0
        lim.xmax = group.attrs.get("xmax")/1000.0
        lim.qmin = ffset.attrs.get("qmin")
        lim.qmax = ffset.attrs.get("qmax")
        if ( "ymin" in group.attrs ):
            lim.ymin = group.attrs.get("ymin")/1000.0
        else:
            lim.ymin = lim.xmin
        if ( "ymax" in group.attrs ):
            lim.ymax = group.attrs.get("ymax")/1000.0
        else:
            lim.ymax = lim.xmax

        ffPlot.uid = group.attrs.get("uid")
        ampPlot.uid = ffPlot.uid
        phasePlot.uid = ffPlot.uid
        intensityPlot.uid = ffPlot.uid

    ampPlot.limits = lim
    ffPlot.limits = lim
    phasePlot.limits = lim
    intensityPlot.limits = lim
    #ffPlot.cmap = "bone"

    # Compute center of mass
    ffPlot.centerOfMass()
    phasePlot.centerOfMass()
    intensityPlot.centerOfMass()

    root = tk.Tk()
    control = cg.Control( root )
    control = ffPlot.projectionXY( control )
    control = ffPlot.projectionX( control )
    control = ffPlot.projectionY( control )

    control = phasePlot.projectionXY( control )

    control = intensityPlot.projectionXY( control )
    root.mainloop()

if __name__ == "__main__":
    main( sys.argv[1:] )
