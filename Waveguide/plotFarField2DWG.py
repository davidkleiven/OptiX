import sys
from PLOD import controlGUI as cg
import tkinter as tk
sys.path.append("../FresnelFDTD")
sys.path.append("../")
sys.path.append("Scripting/")
import numpy as np
import json
import h5py as h5
from matplotlib import pyplot as plt
import colorScheme as cs
import subprocess
import exactFarFields as eff
DELTA = 4.14E-5 # Salditt et al
BETA = 3.45E-6 # Salditt et al

def main( argv ):
    fname = ""
    MSG = "Usage: python plotFarField2DWG.py --file=<h5file> [--help]\n"
    MSG += "help: Prin this message\n"
    MSG += "file: HDF5 field created by the far field computer.\n"
    MSG += "      Contains one dataset for the exit field and one for the far field\n"
    MSG += "fullExitField: Plot the full exit field as outputted to the HDF5 filed\n"
    MSG += "angles: center,width"
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        elif ( arg.find("--help") != -1 ):
            print (MSG)
            return
        elif ( arg.find("--fullExit") != -1 ):
            fullExitField = True
        elif ( arg.find("--reference=") != -1 ):
            reference = arg.split("--reference=")[1]
        else:
            print ("Unknown argument %s"%(arg))
            return

    root = tk.Tk()
    control = cg.Control( root )
    if ( fname == "" ):
        print ("No file name specified")
        return

    yloc = plt.MaxNLocator(5)
    with h5.File(fname, 'r') as hf:
        group = hf.get("/data")
        exitFieldDset = hf.get("/data/exitField")
        exitField = np.array( exitFieldDset )
        farFieldDset = hf.get("/data/farField")
        wavenumber = group.attrs.get("wavenumber")
        phiMin = farFieldDset.attrs.get("phiMin")
        phiMax = farFieldDset.attrs.get("phiMax")
        farField = np.array( farFieldDset )
        uid = group.attrs.get("uid")
        xmin = group.attrs.get("xmin")
        xmax = group.attrs.get("xmax")

    angle = np.linspace( phiMin, phiMax, len(farField) )

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    # Initialize reference
    color = "#e41a1c"

    ax = fig.add_subplot(1,1,1)
    ax.plot(angle, farField**2, color=color, label="WG")
    ax.set_xlabel("Exit angle (deg)")
    ax.set_ylabel("Intensity (a.u.)")
    ax.set_ylim(bottom=1E-6*np.max(farField**2))

    fname = "Figures/farField%d.svg"%(uid)
    control.attach( fig, ax, fname )

    # Plot the exit field
    x = np.linspace(xmin, xmax, len(exitField))

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(1,1,1)
    ax2 = fig2.add_subplot(1,1,1)
    ax2.plot( x, exitField, color="black")
    ax2.set_xlabel("$x$ (nm)")
    ax2.set_ylabel("Field (a.u.)")

    fname = "Figures/exitField%d.svg"%(uid)
    control.attach( fig2, ax2, fname )
    root.mainloop()

if __name__ == "__main__":
    main( sys.argv[1:])
