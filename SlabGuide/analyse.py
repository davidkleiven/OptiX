USE_LATEX = False
import sys
if USE_LATEX:
    sys.path.append("../FresnelFDTD/")
    import mplLaTeX
    import matplotlib as mpl
    mpl.rcParams.update(mplLaTeX.params)
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
import json

def computeKy( data, fdir ):
      
    y = np.array( data["monitor"]["inside"]["pos"] )
    amp = np.array( data["monitor"]["inside"]["data"] )
    amp /= np.max(amp)
    # Shift y such that center is zero
    y -= y[len(y)/2]
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot( y, amp )
    if ( not USE_LATEX ):
        fig2.show()

    fit = np.arccos(amp)
    fit[y<0.0] = -fit[y<0.0]
    slope, interscept, pvalue, rvalue, stderr = stats.linregress( y, fit )
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot( y, fit, 'ko', fillstyle="none" )
    ax.plot( y, slope*y + interscept, 'k--' )
    if ( USE_LATEX ):
        fig.savefig( fdir+"/kyFit.pdf", bbox_inches="tight" )
    else:
        fig.show()

    print ("Transverse wave vector: %.3E"%(slope))


def plotQFactor( data, fdir ):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    freq = np.array( data["monitor"]["freq"] )
    Q = np.array( data["monitor"]["qFactor"] )
    ax.plot( freq, Q, 'k')
    ax.set_xlabel("Frequency")
    ax.set_ylabel("Quality factor Im($\omega$)/Re($\omega$)")
    if ( USE_LATEX ):
        fname = fdir +"/Qfactor.pdf"
        fig.savefig( fname, bbox_inches="tight")
    else:
        fig.show()
    pos = np.argmax(Q)
    print ("Highest mode frequency: %.4E"%(freq[pos]))
     
def main(argv):
    if ( len(argv) != 2 ):
        print ("Usage: python analyse.py --dfile=<data.json> --fdir=<figdir>")
        return 1

    # Parse arguments
    dfile = ""
    fdir = ""
    for arg in argv:
        if ( arg.find("--dfile=") != -1 ):
            dfile = arg.split("--dfile=")[1]
        elif ( arg.find("--fdir=") != -1 ):
            fdir = arg.split("--fdir=")[1] 
    # Consistency check
    if ( dfile == "" ):
        print ("No data file specified...")
        return 1
    elif ( fdir == "" ):
        print ("No figure directory specified...")
        return 1

    # Read data
    try:
        infile = open( dfile, 'r' )
        data = json.load(infile)
    except:
        print ("Could not open file %s"%(dfile))
        return 1
    # Compute the transverse k-vector
    computeKy( data, fdir )
    plotQFactor( data, fdir )

    if ( not USE_LATEX ):
        plt.show()

if __name__  == "__main__":
    main(sys.argv[1:])
  
