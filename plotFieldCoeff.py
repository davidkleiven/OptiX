import numpy as np
import mplLaTeX
import matplotlib as mpl
mpl.rcParams.update(mplLaTeX.params)
from matplotlib import pyplot as plt
import json

#FOLDERS_s = ["dataPlane/MultInc20/WithEps", "dataPlane/MultInc45/WithEps", "dataPlane/MultInc75/WithEps", "dataPlane/MultInc85/WithEps"]
FOLDERS_s = ["dataPlane/MultInc20/WithEps"]
FOLDERS_p = ["dataPlane/MultInc20p/WithEps", "dataPlane/MultInc45p/WithEps", "dataPlane/MultInc75p/WithEps", \
"dataPlane/MultInc85p/WithEps"]

def transmittedTheta( theta, n1, n2 ):
    return np.arccos( np.sqrt( 1.0 - (n1*np.sin(theta)/n2)**2 ))

def rs(theta, n1, n2):
    theta = theta*np.pi/180.0
    thetaT = transmittedTheta(theta, n1, n2)
    return ( n1*np.cos(theta) - n2*np.cos(thetaT) )/( n1*np.cos(theta) + n2*np.cos(thetaT) )
    
def main():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    theta = np.linspace(0.0, 90.0, 101)
    n1 = 1.0
    n2 = 1.5
    ax.plot( theta, np.abs(rs(theta, n1, n2)), color='black' ) 
    ax2 = ax.twinx()

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

        reflectionNorm = np.array( data["reflection"]["norm"] )
        reflectionPhase = np.array( data["reflection"]["phase"] )
        
        ax.plot( data["angle"], reflectionNorm, 'o', color='black', ms=msize)
        ax2.plot( data["angle"], reflectionPhase, '.', color='grey', ms=msize)

    ax.set_ylim(-1.0, 1.0)
    ax2.set_ylim(-np.pi/2.0,np.pi/2.0)
    fig.savefig( "Figures/fieldCoefficients.pdf", bbox_inches="tight" )

if __name__ == '__main__':
    main()
