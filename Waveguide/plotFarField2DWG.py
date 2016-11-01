import sys
sys.path.append("../FresnelFDTD")
sys.path.append("../")
import mplLaTeX as ml
import matplotlib as mpl
mpl.rcParams.update(ml.params)
import numpy as np
import json
import h5py as h5
from matplotlib import pyplot as plt
import colorScheme as cs
DELTA = 4.14E-5 # Salditt et al
BETA = 3.45E-6 # Salditt et al

def main( argv ):
    fname = ""
    MSG = "Usage: python plotFarField2DWG.py --file=<h5file> [--help]\n"
    MSG += "help: Prin this message\n"
    MSG += "file: HDF5 field created by the far field computer.\n"
    MSG += "      Contains one dataset for the exit field and one for the far field\n"
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        elif ( arg.find("--help") != -1 ):
            print MSG
            return
        else:
            print ("Unknown argument %s"%(arg))
            return

    if ( fname == "" ):
        print ("No file name specified")
        return

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

    maxpos = np.argmax(farField**2)
    angleCenter = angle[maxpos]
    start = np.argmin( np.abs( angle -angleCenter +dAngle) )
    end = np.argmin( np.abs( angle-angleCenter -dAngle) )
    farField = farField[start:end]
    angle = angle[start:end]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(angle, farField**2, color="black")
    #ax.plot(freq, np.abs(pfft)**2, color="red")
    ax.set_yscale("log")
    #ax.set_xlabel("q (nm$^{-1}$)")
    ax.set_xlabel("Exit angle (deg)")
    ax.set_ylabel("Intensity (a.u.)")
    fname = "Figures/farField%d.pdf"%(uid)
    fig.savefig(fname, bbox_inches="tight")
    print ("Figure written to %s"%(fname))

    # Plot the exit field
    indx = np.argmax(np.abs(exitField))
    x = np.linspace(xmin, xmax, len(exitField))
    delta = int( 0.3*len(exitField) )
    minIndx = indx-delta
    maxIndx = indx+delta
    if ( minIndx < 0 ):
        minIndx = 0
    if ( maxIndx >= len(exitField)):
        maxIndx = -1
    exitField = exitField[minIndx:maxIndx]
    x = x[minIndx:maxIndx]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot( x, exitField, color="black")
    ax.set_xlabel("$x$ (nm)")
    ax.set_ylabel("Field (a.u.)")
    fname = "Figures/exitField%d.jpeg"%(uid)
    fig.savefig(fname, bbox_inches="tight")
    print ("Figure written to %s"%(fname))

if __name__ == "__main__":
    main( sys.argv[1:])
