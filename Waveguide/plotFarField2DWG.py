import sys
sys.path.append("../FresnelFDTD")
sys.path.append("../")
#import mplLaTeX as ml
import matplotlib as mpl
#mpl.rcParams.update(ml.params)
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 28
import numpy as np
import json
import h5py as h5
from matplotlib import pyplot as plt
import colorScheme as cs
import subprocess
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
    fullExitField = False
    anglesGiven = False
    angCenter = 0.0
    angWidth = 0.0
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

    angMin = -100.0
    angMax = 100.0
    happy = False
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    print ("Type any non numeric entry to quit")
    while ( not happy ):
        fig.clf()
        ax = fig.add_subplot(1,1,1)
        start = np.argmin( np.abs(angle-angMin))
        end = np.argmin( np.abs(angle-angMax))
        ff = farField[start:end]
        ang = angle[start:end]
        ax.plot(ang, ff**2, color="black")
        ax.set_yscale("log")
        ax.set_xlabel("Exit angle (deg)")
        ax.set_ylabel("Intensity (a.u.)")
        ax.set_ylim(bottom=1E-6*np.max(ff**2))
        plt.show( block=False )
        angMin = raw_input("Min angle:")
        angMax = raw_input("Max angle:")
        try:
            angMin = float(angMin)
            angMax = float(angMax)
        except:
            happy = True
    fname = "Figures/farField%d.svg"%(uid)
    psname = "Figures/farField%d.ps"%(uid)
    fig.savefig(fname, bbox_inches="tight")
    subprocess.call(["inkscape", "--export-ps=%s"%(psname), "--export-latex", fname])
    print ("Figure written to %s"%(fname))
    plt.close()

    # Plot the exit field
    indx = np.argmax(np.abs(exitField))
    x = np.linspace(xmin, xmax, len(exitField))

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    happy = False
    start = 0
    end = -1
    while ( not happy ):
        fig.clf()
        ax = fig.add_subplot(1,1,1)
        xx = x[start:end]
        ef = exitField[start:end]
        ax.plot( xx, ef, color="black")
        ax.set_xlabel("$x$ (nm)")
        ax.set_ylabel("Field (a.u.)")
        plt.show( block=False )
        xmin = raw_input("Minval: ")
        xmax = raw_input("Maxval: ")
        try:
            xmin = float(xmin)
            xmax = float(xmax)
            start = np.argmin( np.abs(x-xmin))
            end = np.argmin( np.abs(x-xmax))
        except:
            happy = True

    fname = "Figures/exitField%d.svg"%(uid)
    fig.savefig(fname, bbox_inches="tight")
    psname = "Figures/exitField%d.ps"%(uid)
    subprocess.call(["inkscape", "--export-ps=%s"%(psname), "--export-latex", fname])
    print ("Figure written to %s"%(fname))

if __name__ == "__main__":
    main( sys.argv[1:])
