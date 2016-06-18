import numpy as np
import mplLaTeX
import matplotlib as mpl
mpl.rcParams.update(mplLaTeX.params)
from matplotlib import pyplot as plt
import json

FOLDERS_s = ["dataPlane/MultInc5/WithEps", "dataPlane/MultInc20/WithEps", "dataPlane/MultInc45/WithEps", \
"dataPlane/MultInc75/WithEps", "dataPlane/MultInc85/WithEps"]
#FOLDERS_s = ["dataPlane/MultInc20/WithEps"]
FOLDERS_p = ["dataPlane/MultInc5p/WithEps", "dataPlane/MultInc20p/WithEps", "dataPlane/MultInc45p/WithEps", \
"dataPlane/MultInc75p/WithEps", "dataPlane/MultInc85p/WithEps"]
#FOLDERS_p = ["dataPlane/MultInc20p/WithEps"]
def transmittedTheta( theta, n1, n2 ):
    return np.arccos( np.sqrt( 1.0 - (n1*np.sin(theta)/n2)**2 ))

def rs(theta, n1, n2):
    theta = theta*np.pi/180.0
    thetaT = transmittedTheta(theta, n1, n2)
    return ( n1*np.cos(theta) - n2*np.cos(thetaT) )/( n1*np.cos(theta) + n2*np.cos(thetaT) )

def rp(theta, n1, n2): 
    theta = theta*np.pi/180.0
    thetaT = transmittedTheta(theta, n1, n2)
    return ( n1*np.cos(thetaT) - n2*np.cos(theta) )/( n1*np.cos(thetaT) + n2*np.cos(theta) )
    
def main():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    figPhase = plt.figure()
    axPhase = figPhase.add_subplot(111)
    theta = np.linspace(0.0, 90.0, 101)
    n1 = 1.0
    n2 = 1.5
    ax.plot( theta, np.abs(rs(theta, n1, n2)), color='black' ) 
    ax.plot( theta, np.abs(rp(theta, n1, n2)), color='black' ) 

    msize = 2
    step = 5
    for i in range(0, len(FOLDERS_s)):
        fname_s = FOLDERS_s[i]+'/realFieldFourier.json'
        fname_p = FOLDERS_p[i]+'/realFieldFourier.json'
        try:
            infile = open( fname_s, 'r')
            data = json.load(infile)
            infile.close()
        except:
            print ("Could not open file %s"%(fname_s))
            continue  

        try:
            infile = open( fname_p, 'r' )
            data_p = json.load(infile)
            infile.close()
        except:
            print ("Could not open file %s"%(fname_p))
            continue

        reflectionNorm = np.array( data["reflection"]["norm"] )[0:-1:step]
        reflectionPhase = np.array( data["reflection"]["phase"] )[0:-1:step]
        angle = np.array( data["angle"] )[0:-1:step]

        refPNorm = np.array( data_p["reflection"]["norm"] )[0:-1:step]
        refPPhase = np.array( data_p["reflection"]["phase"])[0:-1:step]
        angle_p = np.array( data_p["angle"] )[0:-1:step]
        
        ax.plot( angle, reflectionNorm, 'o', color='black', ms=msize)
        ax.plot( angle_p, refPNorm, '^', color='black', ms=msize)
        axPhase.plot( angle, reflectionPhase, '.', color='black', ms=msize)

    ax.set_ylim(0.0, 1.0)
    axPhase.set_ylim(-np.pi,np.pi)

    # Set proper labels on the phase plot
    ymarks = np.array([-np.pi, -np.pi/2.0, 0.0, np.pi/2.0, np.pi])
    axPhase.set_yticks(ymarks)
    ytics = ["$-\pi$", "$\\frac{\pi}{2}$", "0", "$\\frac{\pi}{2}$", "$\pi$"]
    axPhase.set_yticklabels(ytics)
    fig.savefig( "Figures/fieldCoefficients.pdf", bbox_inches="tight" )
    figPhase.savefig( "Figures/fieldCoefficientPhase.pdf", bbox_inches="tight" )

if __name__ == '__main__':
    main()
