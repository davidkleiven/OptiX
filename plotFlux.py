import numpy as np
import mplLaTeX
import matplotlib as mpl
mpl.rcParams.update(mplLaTeX.params)
from matplotlib import pyplot as plt
import json
import sys

FOLDERS_s = ["dataPlane/MultInc20s/WithEps", "dataPlane/MultInc45s/WithEps", "dataPlane/MultInc75s/WithEps", \
"dataPlane/MultInc85s/WithEps"]
FOLDERS_p = ["dataPlane/MultInc20p/WithEps", "dataPlane/MultInc45p/WithEps", "dataPlane/MultInc75p/WithEps", \
"dataPlane/MultInc85p/WithEps"]
SUBDIR="WithEps"
POLARISATIONS=["s", "p"]

def transmissionAngle(theta, n1, n2):
    return np.arccos(np.sqrt( 1.0 - (n1*np.sin(theta)/n2)**2))

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
        print "Usage: python plotFlux.py <datadirbase> <figuredir> <angles>"
        print "incident angle and polarisation and subdir %s will be appended"%(SUBDIR)
        print "Examples: datadirbase = dataplane/MultInc ---> "
        print "s polarisation with 20 deg incidence is located int dataplane/MulcInc20s/%s"%(SUBDIR)
    fig = plt.figure()
    ax = fig.add_subplot(111) 
    figError = plt.figure()
    axError = figError.add_subplot(111)

    theta = np.linspace(0.0, 90.0, 101)
    n1 = 1.0
    n2 = 1.5
    ax.plot(theta, Rs(theta, n1, n2), color='black')
    ax.plot(theta, Rp(theta, n1, n2), color='black')
    ax.plot(theta, Ts(theta, n1, n2), color='black')
    ax.plot(theta, Tp(theta, n1, n2), color='black')
    ax.set_xlabel("Incident angle (deg)")
    ax.set_ylabel("Transmitance/Reflectance")
    folderBase = argv[0]
    fdir = argv[1]
    try:
        incidentAngles = np.array( argv[2:] ).astype(np.int32)
    except:
        print "[PlotFlux] Error when converting arguments to int"
        return

    # Plot simulation results
    step = 10
    hasLabel = False
    markersize=2
    for pol in POLARISATIONS:
        for angle in incidentAngles:
            folder = folderBase+"%d%s/%s"%(angle, pol, SUBDIR)
            try:
                infile = open(folder+"/transmittedFluxNorm.json", 'r')
            except:
                print ("Could not open %s"%(folder+"/transmittedFluxNorm.json"))
                continue
            data = json.load(infile)
            infile.close()

            angles = np.array( data["incidentAngle"] )[0:-1:step]
            T = np.array( data["transmitted"] )[0:-1:step]
            R = np.array( data["reflected"] )[0:-1:step]
            errorR = np.abs(R - Rs(angles, n1, n2))/Rs(angles, n1, n2)
            errorT = np.abs(T - Ts(angles, n1, n2))/Ts(angles, n1, n2)
            if ( hasLabel ):
                ax.plot(angles, R, '^', color='black', ms=markersize)
                ax.plot(angles, T,  'o', color='black', ms=markersize)
                axError.plot(angles, errorR, '^', color='black', ms=markersize)
                axError.plot(angles, errorT, 'o', color='black', ms=markersize)
            else:
                ax.plot(angles, R, '^', color='black', label="$R_\mathrm{%s}$"%(pol), ms=markersize)
                ax.plot(angles, T,  'o', color='black', label="$T_\mathrm{%s}$"%(pol), ms=markersize)
                axError.plot(angles, errorR, '^', color='black', label="$R_\mathrm{%s}$"%(pol), ms=markersize)
                axError.plot(angles, errorT,  'o', color='black', label="$T_\mathrm{%s}$"%(pol), ms=markersize)
                hasLabel = True    
        hasLabel = False
    ax.legend(loc='center left', frameon=False)
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
