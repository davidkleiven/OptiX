from PLOD import controlGUI as cg
import multiWGPlot as mp
import h5py as h5
import sys
import tkinter as tk
import numpy as np
sys.path.append("Scripting")

def main( argv ):
    fname = ""
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        elif ( arg.find("--help") != -1 ):
            print ("Usage: python plotMultiWG.py --file=<h5file> [--help]")
            print ("help: Print this message")
            print ("file: HDF5 file generated by the C++ application multiWG.out")
            return 0
        else:
            print ("Unknown argument %s"%(arg))
            return 1

    if ( fname == "" ):
        print ("No file given!")
        return 1

    multiplot = mp.MultiWGPlot()
    with h5.File( fname, 'r' ) as hf:
        group = hf.get("/data")
        multiplot.zmin = group.attrs.get("zmin")
        multiplot.zmax = group.attrs.get("zmax")
        multiplot.xmin = group.attrs.get("xmin")
        multiplot.xmax = group.attrs.get("xmax")
        multiplot.R = np.array( group.attrs.get("radius") )
        multiplot.angles = np.array( group.attrs.get("angles"))
        uid = group.attrs.get("uid")

        intensity = np.array( hf.get("/data/amplitude") )
        transmittivity = np.array( hf.get("/data/transmittivity") )
        if ( "image" in group.attrs.keys() ):
            multiplot.miniatyrGeo = group.attrs.get( "image" )

    multiplot.miniatyr = False
    # Plot the datasets
    root = tk.Tk()
    ctl = cg.Control( root )
    fig, ax = multiplot.plotIntensity( intensity )
    ctl.attach( fig, ax, "Figures/multiWGIntensity%d.svg"%(uid))

    figLog, axLog = multiplot.plotIntensity( intensity, scale="log" )
    ctl.attach( figLog, axLog, "Figures/multiWGIntensityLog%d.svg"%(uid))

    figT, axT = multiplot.plotTransmittivity( transmittivity )
    ctl.attach( figT, axT, "Figures/multiWGTransmittivity%d.svg"%(uid))
    root.mainloop()

if __name__ == "__main__":
    main( sys.argv[1:] )
