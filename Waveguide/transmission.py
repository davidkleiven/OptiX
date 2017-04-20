import sys
sys.path.append("../FresnelFDTD")
sys.path.append("../")
import mplLaTeX as ml
import matplotlib as mpl
#mpl.rcParams.update(ml.params)
import numpy as np
import h5py as h5
import json
from matplotlib import pyplot as plt
from scipy import stats
import colorScheme as cs

class Dataset:
    def __init__(self):
        self.z = None
        self.data = None
        self.label = ""

class ComparePlots:
    def __init__(self):
        self.datasets = []
        self.figname = "defaultfigure.pdf"

    def addData( self, z, trans, label ):
        dset = Dataset()
        dset.z = z
        dset.data = trans
        dset.label = label
        self.datasets.append(dset)

class Transmission(ComparePlots):
    def __init__(self):
        ComparePlots.__init__(self)

    def plot(self):
        if ( len(self.datasets) == 0 ):
            print ("No datasets.")
            return

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        mxColor = len(cs.COLORS)
        indx = 1
        for dset in self.datasets:
            ax.plot( dset.z/1E6, np.log(dset.data), label=dset.label, color=cs.COLORS[indx%mxColor])
            indx += 1
        ax.set_xlabel("$z$ (mm))")
        ax.set_ylabel("$\ln T$")
        ax.legend(loc="upper right", frameon=False)
        fig.savefig(self.figname, bbox_inches="tight")
        print ("Figure written to %s"%(self.figname))

class ExitFields(ComparePlots):
    def __init__(self):
        ComparePlots.__init__(self)
        self.widthFraction = 0.3

    def cutData( self ):
        for dset in self.datasets:
            maxPos = np.argmax( np.abs(dset.data) )
            N = len( dset.data )
            upperLimit = int(maxPos+self.widthFraction*N)
            lowerLimit = int(maxPos-self.widthFraction*N)
            if ( lowerLimit < 0 ):
                lowerLimit = 0
            if ( upperLimit > N ):
                upperLimit = -1
            dset.data = dset.data[lowerLimit:upperLimit]
            dset.z = dset.z[lowerLimit:upperLimit]

    def plot( self ):
        if ( len(self.datasets) == 0 ):
            print ("No datasets!")
            return
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        mxColor = len(cs.COLORS)
        indx = 0
        for dset in self.datasets:
            indx += 1
            ax.plot( dset.z, dset.data, label=dset.label, color=cs.COLORS[indx%mxColor])
        ax.set_xlabel("x (nm)")
        ax.set_ylabel("Amplitude (a.u.)")
        ax.legend(loc="upper right", frameon=False)
        fig.savefig(self.figname, bbox_inches="tight")
        print ("Figure written to %s"%(self.figname))

def plotTransmission( data, zmin, zmax, uid, control ):
    z = np.linspace( zmin, zmax, len(data))
    fitStart = int( 3*len(data)/4 )
    slope, interscept, rvalue, pvalue, stderr = stats.linregress( z[fitStart:], np.log(data[fitStart:]) )
    zFit = np.linspace( np.min(z), np.max(z), 101 )
    fit = interscept + slope*zFit
    print ("Decay length: %.2E mm"%(-1.0/(1E6*slope)))
    print ("Interscept: %.2E"%(np.exp(interscept)))
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot( z/1E6, np.log(data), color="#e41a1c", label="Num" )
    ax.plot( zFit/1E6, fit, color="#377eb8", label="Fit")
    ax.set_xlabel("$z$ (mm)")
    ax.set_ylabel("$\ln T$")
    ax.legend(loc="best", frameon=False, labelspacing=0.05)
    fname = "Figures/transmission%d.svg"%(uid)
    control.attach( fig, ax, fname )
