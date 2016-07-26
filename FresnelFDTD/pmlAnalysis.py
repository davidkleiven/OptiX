USE_LATEX = True
import numpy as np
import matplotlib
if ( USE_LATEX ):
    import mplLaTeX as mltx
    matplotlib.rcParams.update(mltx.params)
from matplotlib import pyplot as plt
from scipy import signal
import json
import sys

DDIR = "dataPlane/PMLAnalysis/bkg"
FILES = ["realField5s.json", "realField60s.json", "realField85s.json"]
if ( USE_LATEX ):
    LABELS = ["$\SI{5}{\degree}$", "$\SI{60}{\degree}$", "$\SI{85}{\degree}$", "$\SI{88}{\degree}$"]
else:
    LABELS = ["5 deg", "60 deg", "85 deg", "88 deg"]

MARKER = ["o", "^", "s", "p"]
FDIR = "Figures"

# Add directory
for i in range(0,len(FILES)):
    FILES[i] = DDIR+"/"+FILES[i]

def main():
    assert ( len(FILES) <= len(LABELS) )
    assert ( len(FILES) <= len(MARKER) )

    fig = plt.figure()  
    ax = fig.add_subplot(111)
    for i in range(0, len(FILES)):
        try:
            infile = open(FILES[i], 'r')
        except:
            print ("Could not open file %s"%(FILES[i]))
            continue
        
        data = json.load( infile )
        infile.close()
        amplitude = np.array( data["PMLMonitors"]["MaxAmplitude"] )
        position = np.array( data["PMLMonitors"]["position"] )
        position = position[::-1]
        amplitude /= np.max( amplitude )
        ax.plot( position, amplitude, color="black", label=LABELS[i], marker=MARKER[i], fillstyle="none" )

    ax.set_xlabel("Distance inside PML layer")
    ax.set_ylabel("Amplitude")
    ax.legend( loc="lower left", frameon=False )
    
    if ( USE_LATEX ):
        fname = FDIR +"/pmlDamping.pdf"
        fig.savefig( fname, bbox_inches="tight" )
        print ("Figure saved to %s"%(fname))
    else:
        plt.show()

if __name__ == "__main__":
    main()

        
