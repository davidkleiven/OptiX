import numpy as np
import mplLaTeX
import matplotlib
matplotlib.rcParams.update(mplLaTeX.params)
from matplotlib import pyplot as plt
from scipy import stats
from scipy import optimize
import json
ANGLES = [5, 20, 45, 75, 85]
POLARISATIONS=["s","p"]
BASE = "dataPlane/MultInc"
SUBDIR = "WithEps"
FNAME = "realFieldFourier.json"

def findReflectionAngle(theta_r, theta_in, waveNumber, distanceFromSlab, phase):
    # Force the solver to stay within +- pi
    if (( theta_r > np.pi/2.0 ) or ( theta_r < 0.0 )):
        return np.inf
    return waveNumber*distanceFromSlab*(np.sin(theta_in)*(np.tan(theta_in)+np.tan(theta_r)) \
    - 1.0/np.cos(theta_in) - 1.0/np.cos(theta_r)) - phase

def expectedPhase(freq, fcenter, thetaCenter, distanceFromPlane):
    sqrtArg = freq**2 - (fcenter*np.sin(thetaCenter))**2 
    sqrtArg[sqrtArg < 0.0] = 0.0
    return 2.0*np.pi*2.0*distanceFromPlane*np.sqrt( sqrtArg )

def brewster(n1, n2):
    return np.arctan(n2/n1)

def main():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    tanReflTimesTanAngle = None
    tanAngle = None
    step = 10
    brewsterAngle = -1.0
    for pol in POLARISATIONS:
        hasLabel = False
        for theta in ANGLES:
            fname = BASE+"%d%s/"%(theta, pol)+SUBDIR +"/"+FNAME
            try:
                infile = open(fname,'r')
                data = json.load(infile)
                infile.close()
            except:
                print ("Could not find file %s"%(fname))
                continue

            if ( brewsterAngle < 0.0 ):
                n1 = np.sqrt( data["geometry"]["EpsilonLow"])
                n2 = np.sqrt( data["geometry"]["EpsilonHigh"] )
                brewsterAngle = brewster(n1, n2)
            distanceFromPlane = data["reflected"]["position"] - data["geometry"]["slabPosition"]
            phase = np.array( data["reflected"]["phase"] )
            angle = np.array( data["reflected"]["angle"] )*np.pi/180.0
            k = 2.0*np.pi*np.array( data["reflected"]["frequency"] )

            angle = angle[::step]
            k = k[::step]
            phase = phase[::step]
            x = 2.0*distanceFromPlane*k*np.cos(angle)
            m = np.floor( x/(2.0*np.pi) )
            phase -= m*2.0*np.pi
            marker = '.'
            label="$\\angle r_\mathrm{s}$"
            msize = 5
            if ( pol == "p" ):
                marker = 'x'
                msize=2
                label="$\\angle r_\mathrm{p}$"
            if ( hasLabel ):
                ax.plot(x, phase, marker, color="black", markersize=msize, fillstyle="none")
            else: 
                ax.plot(x, phase, marker, color="black", markersize=msize, fillstyle="none", label=label)
                hasLabel = True

    x1, x2 = ax.get_xlim()
    x = np.linspace(0.9*x1, 1.1*x2, 11)
    ax.plot( x, -x+np.pi, color='black')
    ax.plot( x, -x, color='black')
    ax.set_xlabel("$\\frac{2y\omega}{c} \cos \\theta_\mathrm{i}$")
    ax.set_ylabel("$\\phi_\mathrm{\omega}$")
    ax.legend(loc="upper right", frameon=False)
    fig.savefig("Figures/tanReflection.pdf", bbox_inches="tight")

if __name__ == "__main__":
    main()
