import numpy as np
import mplLaTeX
import matplotlib as mpl
mpl.rcParams.update(mplLaTeX.params)
from matplotlib import pyplot as plt
import json

FOLDERS_s = ["dataPlane/MultInc20/WithEps", "dataPlane/MultInc45/WithEps", "dataPlane/MultInc75/WithEps", "dataPlane/MultInc85/WithEps"]
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

        angles = data["incidentAngle"][0:-1:step]
        T = data["transmitted"][0:-1:step]
        R = data["reflected"][0:-1:step]
        if ( hasLabel ):
            ax.plot(angles, R, '^', color='black', ms=markersize)
            ax.plot(angles, T,  'o', color='black', ms=markersize)
        else:
            ax.plot(angles, R, '^', color='black', label="$R_s$", ms=markersize)
            ax.plot(angles, T,  'o', color='black', label="$T_s$", ms=markersize)
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

        angles = data["incidentAngle"][0:-1:step]
        T = data["transmitted"][0:-1:step]
        R = data["reflected"][0:-1:step]
        if ( hasLabel ):
            ax.plot(angles, R, 's', color='black', ms=markersize)
            ax.plot(angles, T,  'x', color='black', ms=markersize)
        else:
            ax.plot(angles, R, 's', color='black', label="$R_p$", ms=markersize)
            ax.plot(angles, T, 'x', color='black', label="$T_p$", ms=markersize)
            hasLabel=True
    ax.legend(loc='center left', frameon=False)
    fig.savefig("Figures/powerCoefficients.pdf", bbox_inches="tight")

if __name__ == '__main__':
    main() 
