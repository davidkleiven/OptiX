import sys
sys.path.append("../FresnelFDTD")
import mplLaTeX as ml
import matplotlib as mpl
mpl.rcParams.update(ml.params)
import numpy as np
import h5py as h5
import json
from matplotlib import pyplot as plt
from scipy import stats

def plotTransmission( data, stat ):
    z = np.linspace( stat["Transmission"]["zStart"], stat["Transmission"]["zEnd"], len(data))
    fitStart = int( 3*len(data)/4 )
    slope, interscept, rvalue, pvalue, stderr = stats.linregress( z[fitStart:], np.log(data[fitStart:]) )
    zFit = np.linspace( np.min(z), np.max(z), 101 )
    fit = np.exp(interscept)*np.exp(slope*zFit)
    print ("Decay length: %.2E mm"%(-1.0/(1E6*slope)))
    print ("Interscept: %.2E"%(np.exp(interscept)))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot( z/1E6, data, 'k.', ms=4, fillstyle="none" )
    ax.plot ( zFit/1E6, fit, 'k--')
    ax.set_xlabel("$z (\SI{}{\milli\meter}$)")
    ax.set_ylabel("Tranmission")
    ax.set_yscale("log")
    fname = "Figures/transmission%d.pdf"%(stat["UID"])
    fig.savefig(fname, bbox_inches="tight")
    print ("Figure written to %s"%(fname))
