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
    intensityPlot = proj3D.Intensity()

    with h5.File( fname, 'r' ) as hf:
        group = hf.get("/data")
        ffset = hf.get("/data/farField")
        ffPlot.data = np.array( ffset )
        phasePlot.data = np.array( hf.get("/data/exitPhase") )
        intensityPlot.data = np.array( hf.get("/data/exitIntensity") )
        lim.zmin = group.attrs.get("zmin")
        lim.zmax = group.attrs.get("zmax")
        lim.xmin = group.attrs.get("xmin")/1000.0
        lim.xmax = group.attrs.get("xmax")/1000.0
        lim.qmin = ffset.attrs.get("qmin")
        lim.qmax = ffset.attrs.get("qmax")
        lim.ymin = lim.xmin
        lim.ymax = lim.xmax
        ffPlot.uid = group.attrs.get("uid")

    ampPlot.limits = lim
    ffPlot.limits = lim
    phasePlot.limits = lim
    intensityPlot.limits = lim
    ffPlot.cmap = "bone"

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
    control = phasePlot.projectionX( control )
    control =phasePlot.projectionY( control )

    control = intensityPlot.projectionXY( control )
    root.mainloop()

if __name__ == "__main__":
    main( sys.argv[1:] )
