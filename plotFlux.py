import numpy as np
import mplLaTeX
import matplotlib as mpl
mpl.rcParams.update(mplLaTeX.params)
from matplotlib import pyplot as plt
import json

FOLDERS_s = ["dataPlane/MultInc20s/WithEps", "dataPlane/MultInc45s/WithEps", "dataPlane/MultInc75s/WithEps", \
"dataPlane/MultInc85s/WithEps"]
FOLDERS_p = ["dataPlane/MultInc20p/WithEps", "dataPlane/MultInc45p/WithEps", "dataPlane/MultInc75p/WithEps", \
"dataPlane/MultInc85p/WithEps"]

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

def main():
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

    # Plot simulation results
    step = 10
    hasLabel = False
    markersize=2
    for folder in FOLDERS_s:
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
            ax.plot(angles, R, '^', color='black', label="$R_\mathrm{s}$", ms=markersize)
            ax.plot(angles, T,  'o', color='black', label="$T_\mathrm{s}$", ms=markersize)
            axError.plot(angles, errorR, '^', color='black', label="$R_\mathrm{s}$", ms=markersize)
            axError.plot(angles, errorT,  'o', color='black', label="$T_\mathrm{s}$", ms=markersize)
            hasLabel = True
        
    hasLabel = False
    for folder in FOLDERS_p:
        try:
            infile = open(folder+"/transmittedFluxNorm.json", 'r')
        except:
            print ("Could not open %s"%(folder+"/transmittedFluxNorm.json"))
            continue
        data = json.load(infile)
        infile.close()

        angles = np.array(data["incidentAngle"])[0:-1:step]
        T = np.array( data["transmitted"] )[0:-1:step]
        R = np.array( data["reflected"] )[0:-1:step]
        errorR = np.abs(R - Rp(angles, n1, n2))/Rp(angles, n1, n2)
        errorT = np.abs(T - Tp(angles, n1, n2))/Tp(angles, n1, n2)
        if ( hasLabel ):
            ax.plot(angles, R, 's', color='black', ms=markersize)
            ax.plot(angles, T,  'x', color='black', ms=markersize)
            axError.plot(angles, errorR, 's', color='black', ms=markersize)
            axError.plot(angles, errorT, 'x', color='black', ms=markersize)
        else:
            ax.plot(angles, R, 's', color='black', label="$R_\mathrm{p}$", ms=markersize)
            ax.plot(angles, T, 'x', color='black', label="$T_\mathrm{p}$", ms=markersize)
            axError.plot(angles, errorR, 's', color='black', label="$R_\mathrm{p}$", ms=markersize)
            axError.plot(angles, errorT,  'x', color='black', label="$T_\mathrm{p}$", ms=markersize)
            hasLabel=True
    ax.legend(loc='center left', frameon=False)
    axError.legend( loc='upper left', frameon=False)
    axError.set_xlabel("Incident angle (deg)")
    axError.set_ylabel("Relative error")
    axError.set_yscale('log')
    fig.savefig("Figures/powerCoefficients.pdf", bbox_inches="tight")
    figError.savefig("Figures/powerCoefficientsError.pdf", bbox_inches="tight")

if __name__ == '__main__':
    main() 
