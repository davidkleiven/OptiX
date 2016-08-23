import sys
sys.path.append("../FresnelFDTD/")
import mplLaTeX
import matplotlib as mpl
mpl.rcParams.update(mplLaTeX.params)
from matplotlib import pyplot as plt
import json
import numpy as np
import fresnelExact as fe

def expectedAttenuationLength( k, eps1, eps2, incangle ):
    n1 = np.sqrt(eps1)
    n2 = np.sqrt(eps2)
    kxi = k*np.sin( incangle*np.pi/180.0 )
    kzt = np.sqrt( k**2 - (n1*kxi/n2)**2 + 0j )
    d = 1.0/(n2.imag*kzt.real + n2.real*kzt.imag)
    return d

# ==== PARAMETERS USED IF NICE PLOT IS TRUE ====
MAX_NUMBER_OF_LINES = 4
COLORS = ["#ca0020", "#f4a582", "#92c5de", "#0571b0"]

def main(argv):
    if ( len(argv) < 1 or len(argv) > 2):
        print ("Usage: python absorption.py --data=<datafilename> [--nice]")
        print ("Option nice: Creates a nice plot with 4 lines and nice colors")
        return 1

    try:
        fname = argv[0].split("--data=")[1]
    except:
        print ("Error when parsing command line argument...")
        return 1

    with open(fname, 'r') as infile:
        data = json.load(infile)

    niceplot = False
    for arg in argv:
        if (arg.find("--nice")):
            niceplot = True
    x = np.array( data["absorption"]["position"] )
    values = data["absorption"]["absMonitor"]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    smallestvalue = 1E16
    if (niceplot):
        delta = int(len(values)/(2.0*(MAX_NUMBER_OF_LINES-1)))
    else:
        delta = 1
    counter = 0
    nextIndx = 0
    plotNumber = 0
    for entry in values:
        if ( entry["polarisation"] == "p" ):
            continue
        if ( counter<nextIndx ):
            counter += 1
            continue
        k = 0.01
        eps1 = 1.0
        eps2 = 0.99998+1.99998E-6*1j
        angle = entry["angle"]
        attenuation = expectedAttenuationLength( k, eps1, eps2, angle )
        y = np.abs( fe.ts(eps1, eps2, 1.0, 1.0, angle, k)*np.exp(-x/attenuation) )**2
        amp = np.array(entry["amplitude"])
        wavelength = 2.0*np.pi/k
        if ( niceplot ):
            ax.plot( x/wavelength, amp, 'o', fillstyle="none", color="black", ms=3)
            ax.plot( x/wavelength, y, color=COLORS[plotNumber], label="$\SI{%2.2f}{\degree}$"%(angle))
            xNorm = x/wavelength
            plotNumber += 1
        else:
            ax.plot( x/wavelength, amp, ".")
            ax.plot( x/wavelength,y)
        if (( np.min(amp) < smallestvalue ) and (np.min(entry) > 0.0)):
            smallestvalue = np.min(amp)
        counter += 1
        nextIndx += delta
    ax.set_yscale("log")
    ax.set_ylim(0.8*smallestvalue,10)
    ax.legend( loc="lower left", frameon=False)
    ax.set_xlabel("Distance inside slab ($\lambda$)")
    ax.set_ylabel("TE amplitude ratio")
    fname = "Figures/absorption.pdf"
    fig.savefig( fname, bbox_inches="tight" )
    print ("Figure written to %s"%(fname))

if __name__ == "__main__":
    main(sys.argv[1:])
    
    
