import sys
sys.path.append("../FresnelFDTD/")
import mplLaTeX
import matplotlib as mpl
mpl.rcParams.update(mplLaTeX.params)
from matplotlib import pyplot as plt
import json
import numpy as np

def main(argv):
    if ( len(argv) != 1 ):
        print ("Usage: python absorption.py --data=<datafilename>")
        return 1

    try:
        fname = argv[0].split("--data=")[1]
    except:
        print ("Error when parsing command line argument...")
        return 1

    with open(fname, 'r') as infile:
        data = json.load(infile)

    x = np.array( data["absorption"]["position"] )
    values = data["absorption"]["absMonitor"]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    smallestvalue = 1E16
    for entry in values:
        amp = np.array(entry["amplitude"])
        ax.plot( x, amp, '.')
        if (( np.min(amp) < smallestvalue ) and (np.min(entry) > 0.0)):
            smallestvalue = np.min(amp)
    ax.set_yscale("log")
    ax.set_ylim(0.8*smallestvalue,1)
    fname = "Figures/absorption.pdf"
    fig.savefig( fname, bbox_inches="tight" )

if __name__ == "__main__":
    main(sys.argv[1:])
    
    
