import numpy as np
import mplLaTeX
import matplotlib as mpl
mpl.rcParams.update(mplLaTeX.params)
from matplotlib import pyplot as plt
import json

FOLDERS_s = ["dataPlane/MultInc20/WithEps", "dataPlane/MultInc45/WithEps", "dataPlane/MultInc75/WithEps", "dataPlane/MultInc85/WithEps"]
FOLDERS_p = ["dataPlane/MultInc20p/WithEps", "dataPlane/MultInc45p/WithEps", "dataPlane/MultInc75p/WithEps", \
"dataPlane/MultInc85p/WithEps"]

def main():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax2 = ax.twinx()

    msize = 2
    for folder in FOLDERS_s:
        fname = folder+'/realFieldFourier.json'
        try:
            infile = open( fname, 'r')
        except:
            print ("Could not open file %s"%(fname))
            continue  

        data = json.load(infile)
        infile.close()

        reflectionReal = np.array( data["reflection"]["real"] )
        reflectionImag = np.array( data["reflection"]["imag"] )
        
        reflectionPhase = np.arctan( reflectionImag/reflectionReal )*180.0/np.pi
        reflectionNorm = np.sqrt( reflectionReal**2 + reflectionImag**2 )
        ax.plot( data["angle"], reflectionNorm, 'o', color='black', ms=msize)
        #ax2.plot( data["angle"], reflectionPhase, '.', color='grey')

    fig.savefig( "Figures/fieldCoefficients.pdf", bbox_inches="tight" )

if __name__ == '__main__':
    main()
