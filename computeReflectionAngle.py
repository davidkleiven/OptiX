import numpy as np
import mplLaTeX
import matplotlib
matplotlib.rcParams.update(mplLaTeX.params)
from matplotlib import pyplot as plt
from scipy import stats
from scipy import optimize
import json

FOLDERS = ["dataPlane/MultInc5s/WithEps", "dataPlane/MultInc20s/WithEps", "dataPlane/MultInc45s/WithEps", \
"dataPlane/MultInc75s/WithEps", "dataPlane/MultInc85s/WithEps"]#, "dataPlane/MultInc5p/WithEps", \
#"dataPlane/MultInc20p/WithEps", "dataPlane/MultInc45p/WithEps", \
#"dataPlane/MultInc75p/WithEps", "dataPlane/MultInc85p/WithEps"]

def findReflectionAngle(theta_r, theta_in, waveNumber, distanceFromSlab, phase):
    # Force the solver to stay within +- pi
    if (( theta_r > np.pi/2.0 ) or ( theta_r < 0.0 )):
        return np.inf
    return waveNumber*distanceFromSlab*(np.sin(theta_in)*(np.tan(theta_in)+np.tan(theta_r)) \
    - 1.0/np.cos(theta_in) - 1.0/np.cos(theta_r)) - phase

def main():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    tanReflTimesTanAngle = None
    tanAngle = None
    step = 10
    hasLabel = False
    for folder in FOLDERS:
        fname = folder+"/realFieldFourier.json"
        try:
            infile = open(fname,'r')
            data = json.load(infile)
            infile.close()
        except:
            print ("Could not find file %s"%(fname))
            continue
        distanceFromPlane = data["reflected"]["position"] - data["geometry"]["slabPosition"]
        phase = np.array( data["reflected"]["phase"] )
        angle = np.array( data["reflected"]["angle"] )*np.pi/180.0
        k = 2.0*np.pi*np.array( data["reflected"]["frequency"] )

        angle = angle[::step]
        k = k[::step]
        phase = phase[::step]
        phase[phase<0.0] += 2.0*np.pi
        phase = 2.0*np.pi - phase
        phase /= (2.0*k*distanceFromPlane) 
        if ( hasLabel ):
            ax.plot(np.cos(angle), phase, '.', color="black", markersize=5, fillstyle="none")
        else: 
            ax.plot(np.cos(angle), phase, '.', color="black", markersize=5, fillstyle="none", label="Numerical")
            hasLabel = True
    cosAngle = np.linspace(0.0, 1.0, 11)
    ax.plot(cosAngle, cosAngle, color="black", label="Analytical")
    ax.set_xlabel("$\cos \\theta_\mathrm{i}$")
    ax.set_ylabel("$\\frac{\\delta}{2ky}$", rotation=0)
    ax.legend(loc="upper left", frameon=False)
    fig.savefig("Figures/tanReflection.pdf", bbox_inches="tight")

if __name__ == "__main__":
    main()
