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

    ampPlot.name = "intensity"
    phasePlot.name = "phase"

    with h5.File( fname, 'r' ) as hf:
        group = hf.get("/data")
        ampPlot.data = np.array( hf.get("/data/amplitude") )
        #phasePlot.data = np.array( hf.get("/data/phase") )
        lim.zmin = group.attrs.get("zmin")
        lim.zmax = group.attrs.get("zmax")
        lim.xmin = group.attrs.get("xmin")
        lim.xmax = group.attrs.get("xmax")
        lim.ymin = lim.xmin
        lim.ymax = lim.xmax
        #lim.ymin = group.attrs.get("ymin")
        #lim.ymax = group.attrs.get("ymax")

    ampPlot.limits = lim

    print (ampPlot.data)
    # Compute center of mass
    ampPlot.centerOfMass()
    #phasePlot.centerOfMass()

    root = tk.Tk()
    control = cg.Control( root )
    ampPlot.sliceZ( cg )

    root.mainloop()

if __name__ == "__main__":
    main( sys.argv[1:] )
