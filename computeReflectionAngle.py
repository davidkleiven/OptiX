import numpy as np
import mplLaTeX
import matplotlib
matplotlib.rcParams.update(mplLaTeX.params)
from matplotlib import pyplot as plt
from scipy import stats
from scipy import optimize
import json

FOLDERS = ["dataPlane/MultInc5/WithEps", "dataPlane/MultInc20/WithEps", "dataPlane/MultInc45/WithEps", \
"dataPlane/MultInc75/WithEps", "dataPlane/MultInc85/WithEps", "dataPlane/MultInc5p/WithEps", \
"dataPlane/MultInc20p/WithEps", "dataPlane/MultInc45p/WithEps", \
"dataPlane/MultInc75p/WithEps", "dataPlane/MultInc85p/WithEps"]

def findReflectionAngle(theta_r, theta_in, waveNumber, distanceFromSlab, phase):
    # Force the solver to stay within +- pi
    if ( np.abs(theta_r) > np.pi ):
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
        distanceFromPlane = data["geometry"]["sourcePosition"] - data["geometry"]["slabPosition"]
        phase = np.array( data["reflection"]["phase"] )
        angle = np.array( data["angle"] )*np.pi/180.0
        k = 2.0*np.pi*np.array( data["frequency"] )

        angle = angle[::step]
        k = k[::step]
        phase = phase[::step]
        reflAngle = np.zeros( len(angle) )
        for i in range(0, len(angle)):
            reflAngle[i] = optimize.fsolve( findReflectionAngle, angle[i], args=(angle[i], k[i], distanceFromPlane, phase[i] ))

        if ( hasLabel ):
            ax.plot(angle, reflAngle, '.', color="black", markersize=5, fillstyle="none")
        else: 
            ax.plot(angle, reflAngle, '.', color="black", markersize=5, fillstyle="none", label="Numerical")
            hasLabel = True
    anglePlot = np.linspace(0.0, np.pi/2.0, 11)
    ax.plot(anglePlot, anglePlot, color="black", label="$\\theta_\mathrm{r} = \\theta_\mathrm{i}$")
    ax.set_xlabel("$\\theta_\mathrm{i}$")
    ax.set_ylabel("$\\theta_\mathrm{r}$")
    ax.legend(loc="upper left", frameon=False)
    ticks = [0.0, np.pi/8.0, np.pi/4.0, 3.0*np.pi/8.0, np.pi/2.0]
    labels = ["$0$", "$\\frac{\pi}{8}$", "$\\frac{\pi}{4}$", "$\\frac{3\pi}{8}$", "$\\frac{\pi}{2}$"]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)
    fig.savefig("Figures/tanReflection.pdf", bbox_inches="tight")

if __name__ == "__main__":
    main()
