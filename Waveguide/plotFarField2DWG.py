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
    MSG += "reference: Plot a reference. Options: twoMirrorGaussians, slits"
    fullExitField = False
    anglesGiven = False
    angCenter = 0.0
    angWidth = 0.0
    reference = "default"
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        elif ( arg.find("--help") != -1 ):
            print (MSG)
            return
        elif ( arg.find("--fullExit") != -1 ):
            fullExitField = True
        elif ( arg.find("--angles=") != -1 ):
            anglesGiven = True
            angs = arg.split("--angles=")[1]
            angCenter = float( angs.split(",")[0] )
            angWidth = float( angs.split(",")[1] )
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
        exitFieldDset = hf.get("exitField")
        exitField = np.array( exitFieldDset )
        farFieldDset = hf.get("farField")
        wavenumber = farFieldDset.attrs.get("wavenumber")
        dx = farFieldDset.attrs.get("gridspacing")
        farField = np.array( farFieldDset )
        uid = farFieldDset.attrs.get("uid")
        xmin = exitFieldDset.attrs.get("xmin")
        xmax = exitFieldDset.attrs.get("xmax")

    if ( dx is None ):
        dx = 1.0
    if ( wavenumber is None ):
        wavenumber = 1.0
    if ( uid is None ):
        uid = 0
    if ( xmin is None ):
        xmin = 0.0
    if ( xmax is None ):
        xmax = 100.0

    #pfft = np.fft.fft(exitField)
    #freq = np.fft.fftfreq(len(pfft), d=dx)
    #freq = 2.0*np.pi*np.fft.fftshift(freq)
    fmin = -np.pi
    fmax = np.pi
    q = np.linspace(fmin,fmax,len(farField))/dx
    angle = q/wavenumber
    angle *= 180.0/np.pi
    dAngle = 1.0

    angMin = -100.0
    angMax = 100.0
    happy = False
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    # Initialize reference
    color = "#e41a1c"
    uselog = True
    uselegend = True
    if ( reference == "twoMirrorGaussians" ):
        refPlot = eff.TwoMirrorGaussians()
        refPlot.angle = 0.57*np.pi/180.0
        refPlot.wavenumber = wavenumber
        refPlot.wgwidth = 50.0
    elif ( reference == "slit" ):
        refPlot = eff.YoungSlit()
        refPlot.width = 100.0
        refPlot.separation = 200.0
        uselog = False
    else:
        refPlot = eff.ExactFarFieldDefault()
        color = "black"
        uselegend = False

    ax = fig.add_subplot(1,1,1)
    ax.plot(angle, farField**2, color=color, label="WG")
    refPlot.fit( np.linspace(q[0],q[-1],len(farField)), farField**2 )
    refPlot.normalize( farField**2 )
    ax = refPlot.plot( ax, q[0], q[-1] )
    ax.set_xlabel("Exit angle (deg)")
    ax.set_ylabel("Intensity (a.u.)")
    ax.set_ylim(bottom=1E-6*np.max(farField**2))
    if ( uselegend ):
        ax.legend(loc="upper right", labelspacing=0.01, frameon=False)

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
