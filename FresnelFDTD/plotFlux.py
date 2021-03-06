import numpy as np
import mplLaTeX
import matplotlib as mpl
mpl.rcParams.update(mplLaTeX.params)
from matplotlib import pyplot as plt
import json
import sys

SUBDIR="WithEps"
POLARISATIONS=["s", "p"]

def transmissionAngle(theta, n1, n2):
    sqrtArg = 1.0 - (n1*np.sin(theta)/n2)**2
    sqrtArg[sqrtArg < 0.0] = 0.0
    cosAngle = np.sqrt( sqrtArg )
    return np.arccos( cosAngle ) 

def Rs(theta, n1, n2):
    theta = theta*np.pi/180.0
    thetaT = transmissionAngle(theta, n1, n2)
    return (( n1*np.cos(theta) - n2*np.cos(thetaT) )/( n1*np.cos(theta) + n2*np.cos(thetaT)))**2

def Rp(theta, n1, n2):
    theta = theta*np.pi/180.0
    thetaT = transmissionAngle(theta, n1, n2)
    return ((n1*np.cos(thetaT) - n2*np.cos(theta) )/( n2*np.cos(theta) + n1*np.cos(thetaT) ))**2

def Ts(theta, n1, n2):
    return 1.0 - Rs(theta, n1, n2)

def Tp(theta, n1, n2):
    return 1.0 - Rp(theta, n1, n2)

def main(argv):
    if ( len( argv ) < 2 ):
        print "[plotFlux] Usage: python plotFlux.py <datadirbase> <figuredir> <angles>"
        print "[plotFlux] incident angle and polarisation and subdir %s will be appended"%(SUBDIR)
        print "[plotFlux] Examples: datadirbase = dataplane/MultInc ---> "
        print "[plotFlux] s polarisation with 20 deg incidence is located int dataplane/MulcInc20s/%s"%(SUBDIR)
        return 1
    fig = plt.figure()
    ax = fig.add_subplot(111) 
    figError = plt.figure()
    axError = figError.add_subplot(111)
    havePlottedExact = False

    folderBase = argv[0]
    fdir = argv[1]
    try:
        incidentAngles = np.array( argv[2:] ).astype(np.int32)
    except:
        print "[PlotFlux] Error when converting arguments to int"
        return

    # Plot simulation results
    step = 3
    hasLabel = False
    markersize=2
    for pol in POLARISATIONS:
        if ( pol == "s" ):
            markerR = "^"
            markerT = "o"
            fill = "none"
        else:
            markerR = "x"
            markerT = "."
            fill = "full"
        for angle in incidentAngles:
            folder = folderBase+"%d%s/%s"%(angle, pol, SUBDIR)
            fname = folder + "/transmittedFluxNorm.json"
            try:
                infile = open(fname, 'r')
            except:
                print ("Could not open %s"%(fname))
                continue
            data = json.load(infile)
            infile.close()

            try:
                n1 = np.sqrt( data["geometry"]["EpsilonLow"] )
                n2 = np.sqrt( data["geometry"]["EpsilonHigh"] )
            except:
                print ("Could not find epsilon in file %s using default..."%(fname))
                n1 = 1.0
                n2 = 1.5
            angles = np.array( data["incidentAngle"] )[0:-1:step]
            T = np.array( data["transmitted"] )[0:-1:step]
            R = np.array( data["reflected"] )[0:-1:step]
            if ( pol == 's' ):
                errorR = np.abs(R - Rs(angles, n1, n2))/Rs(angles, n1, n2)
                errorT = np.abs(T - Ts(angles, n1, n2))/Ts(angles, n1, n2)
            elif ( pol == 'p' ):
                errorR = np.abs(R - Rp(angles, n1, n2))/Rp(angles, n1, n2)
                errorT = np.abs(T - Tp(angles, n1, n2))/Tp(angles, n1, n2)

            else:
                print ("Unknown polarisation")
            if ( hasLabel ):
                ax.plot(angles, R, markerR, color='black', ms=markersize)
                ax.plot(angles, T,  markerT, color='black', ms=markersize, fillstyle=fill)
                axError.plot(angles, errorR, markerR, color='black', ms=markersize)
                axError.plot(angles, errorT, markerT, color='black', ms=markersize, fillstyle=fill)
            else:
                ax.plot(angles, R, markerR, color='black', label="$R_\mathrm{%s}$"%(pol), ms=markersize)
                ax.plot(angles, T,  markerT, color='black', label="$T_\mathrm{%s}$"%(pol), ms=markersize, fillstyle=fill)
                axError.plot(angles, errorR, markerR, color='black', label="$R_\mathrm{%s}$"%(pol), ms=markersize)
                axError.plot(angles, errorT,  markerT, color='black', label="$T_\mathrm{%s}$"%(pol), ms=markersize, fillstyle=fill)
                hasLabel = True    
            if ( not havePlottedExact ): 
                theta = np.linspace(0.0, 90.0, 101)
                ax.plot(theta, Rs(theta, n1, n2), color='black')
                ax.plot(theta, Rp(theta, n1, n2), color='black')
                ax.plot(theta, Ts(theta, n1, n2), color='black')
                ax.plot(theta, Tp(theta, n1, n2), color='black')
                ax.set_xlabel("Incident angle (deg)")
                ax.set_ylabel("Transmitance/Reflectance")
                havePlottedExact = True
        hasLabel = False
    ax.legend(loc='center left', frameon=False)
    ax.set_ylim(0.0, 1.0)
    axError.legend( loc='upper left', frameon=False)
    axError.set_xlabel("Incident angle (deg)")
    axError.set_ylabel("Relative error")
    axError.set_yscale('log')
    fnameCoeff = "powerCoefficients.pdf"
    fnameError = "powerCoefficientsError.pdf"
    fig.savefig(fdir+"/"+"%s"%(fnameCoeff), bbox_inches="tight")
    figError.savefig(fdir+"/"+"%s"%(fnameError), bbox_inches="tight")

    print "[PlotFlux] Figure written to %s/%s"%(fdir, fnameCoeff)
    print "[PlotFlux] Figure written to %s/%s"%(fdir, fnameError)

if __name__ == '__main__':
    main(sys.argv[1:]) 
