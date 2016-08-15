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

def computeKy( fname, fdir ):
    try:
        infile = open( fname, 'r' )
        data = json.load(infile)
    except:
        print ("Could not open file %s"%(fname))
        return 
      
    y = np.array( data["monitor"]["inside"]["pos"] )
    # Shift y such that center is zero
    y -= y[len(y)/2]
    amp = np.array( data["monitor"]["inside"]["data"] )
    amp /= np.max(amp)

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

    print ("Transverse wave vector: %.3E"%(1.0/slope))

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

    # Compute the transverse k-vector
    computeKy( dfile, fdir )

    if ( not USE_LATEX ):
        plt.show()

if __name__  == "__main__":
    main(sys.argv[1:])
  
