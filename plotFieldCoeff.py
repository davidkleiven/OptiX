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

def ts(theta, n1, n2):
    return 1.0+rs(theta, n1, n2)

def tp(theta, n1, n2):
    theta = theta*np.pi/180.0
    thetaT = transmittedTheta(theta, n1, n2)
    return 2.0*n1*np.cos(theta)/( n1*np.cos(thetaT) + n2*np.cos(theta) )
    
def main():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    figPhase = plt.figure()
    axPhase = figPhase.add_subplot(111)
    theta = np.linspace(0.0, 90.0, 101)
    n1 = 1.0
    n2 = 1.5
    ax.plot( theta, np.abs(rs(theta, n1, n2))**2, color='black' ) 
    ax.plot( theta, np.abs(rp(theta, n1, n2))**2, color='black' ) 
    ax.plot( theta, ts(theta, n1, n2)**2, color='black' )
    ax.plot( theta, tp(theta, n1, n2)**2, color='black' )

    msize = 2
    step = 10
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

        reflectionNorm = np.array( data["reflected"]["norm"] )[0:-1:step]
        transmissionNorm = np.array( data["transmitted"]["norm"] )[0:-1:step]
        angleR = np.array( data["reflected"]["angle"] )[0:-1:step]
        angleT = np.array( data["transmitted"]["angle"] )[0:-1:step]

        refPNorm = np.array( data_p["reflected"]["norm"] )[0:-1:step]
        transPNorm = np.array( data_p["transmitted"]["norm"] )[0:-1:step]
        angle_pR = np.array( data_p["reflected"]["angle"] )[0:-1:step]
        angle_pT = np.array( data_p["transmitted"]["angle"])[0:-1:step]
        
        if ( i==0 ):
            ax.plot( angleR, reflectionNorm**2, 'o', color='black', ms=msize, fillstyle="none", label="$|r_\mathrm{s}|^2$")
            ax.plot( angle_pR, refPNorm**2, '^', color='black', ms=msize, label="$|r_\mathrm{p}|^2$")
            ax.plot( angleT, transmissionNorm**2, 'x', color='black', ms=msize, label="$|t_\matrhrm{s}|^2$")
            ax.plot( angle_pT, transPNorm**2, '.', color='black', ms=msize, label="$|t_\mathrm{p}|^2$")
        else:
            ax.plot( angleR, reflectionNorm**2, 'o', color='black', ms=msize, fillstyle="none")
            ax.plot( angle_pR, refPNorm**2, '^', color='black', ms=msize)
            ax.plot( angleT, transmissionNorm**2, 'x', color='black', ms=msize)
            ax.plot( angle_pT, transPNorm**2, '.', color='black', ms=msize)

    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("Incident angle (deg)")
    ax.set_ylabel("$|r_\mathrm{s}|^2$, $|r_\mathrm{p}|^2$, $|t_\mathrm{s}|^2$, $|t_\mathrm{p}|^2$")
    ax.legend(loc="upper left", ncol=2, frameon=False)

    # Set proper labels on the phase plot
    fig.savefig( "Figures/fieldCoefficients.pdf", bbox_inches="tight" )

if __name__ == '__main__':
    main()
