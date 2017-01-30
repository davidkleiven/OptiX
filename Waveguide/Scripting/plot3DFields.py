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
    phasePlot = proj3D.Plotter3D()
    ffPlot = proj3D.Plotter3D()

    ampPlot.name = "intensity"
    phasePlot.name = "phase"
    ffPlot.name = "farField"

    with h5.File( fname, 'r' ) as hf:
        group = hf.get("/data")
        ffset = hf.get("/data/farField")
        ffPlot.data = np.array( ffset )
        lim.zmin = group.attrs.get("zmin")
        lim.zmax = group.attrs.get("zmax")
        lim.xmin = group.attrs.get("xmin")
        lim.xmax = group.attrs.get("xmax")
        lim.qmin = ffset.attrs.get("qmin")
        lim.qmax = ffset.attrs.get("qmax")
        lim.ymin = lim.xmin
        lim.ymax = lim.xmax
        ffPlot.uid = group.attrs.get("uid")
        #lim.ymin = group.attrs.get("ymin")
        #lim.ymax = group.attrs.get("ymax")

    ampPlot.limits = lim
    ffPlot.limits = lim
    ffPlot.cmap = "gray"

    # Compute center of mass
    ffPlot.centerOfMass()

    root = tk.Tk()
    control = cg.Control( root )
    control = ffPlot.farField( control )
    control = ffPlot.projectionX( control )
    control = ffPlot.projectionY( control )

    root.mainloop()

if __name__ == "__main__":
    main( sys.argv[1:] )
