import sys
sys.path.append("../../FresnelFDTD")
import numpy as np
import mplLaTeX as ml
import matplotlib as mpl
mpl.rcParams.update(ml.params)
from matplotlib import pyplot as plt
from scipy import interpolate as interp
from scipy import signal as sig
from scipy import stats
import wpdExtract as wpd

class AngleComputer:
    def __init__(self):
        self.wpd = None
        self.smoothedWG = None
        self.wgSlope = None
        self.figFolder = "../Figures"
        self.overViewName = "overview.pdf"
        self.angleFname = "rayAngles.pdf"

    def createWG( self ):
        if ( self.wpd is None ):
            print ("No WebPlotDigitizer object given!")
            return 1

        wg = self.wpd.get("waveguide")
        if ( wg is None ):
            print ("Error when reading the waveguide data!")
            return
        x = np.linspace( np.min(wg.x), np.max(wg.x), len(wg.x) )
        interpolator = interp.interp1d( wg.x, wg.y )
        y = interpolator(x)

        # Smooth data
        wLen = int( len(y)/2 )
        if ( wLen%2 == 0 ):
            wLen += 1
        smoothedY = sig.savgol_filter( y, wLen, 3 )
        self.smoothedWG = wpd.WpdDataset()
        self.smoothedWG.x = x
        self.smoothedWG.y= smoothedY

        deriv = sig.savgol_filter(y, wLen, 3, deriv=1, delta=x[1]-x[0] )
        self.wgSlope = wpd.WpdDataset()
        self.wgSlope.x = x
        self.wgSlope.y = deriv
        self.angles = None
        self.angleX = None

    def plotView( self ):
        if ( self.wpd is None ):
            print ("No WebPlotDigitizer object given!" )
            return 1

        if ( self.smoothedWG is None ):
            self.createWG()

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        for dset in self.wpd.dsets:
            if ( dset.name == "waveguide" ):
                continue
            ax.plot( dset.x, dset.y, color="red", marker='.')

        # Plot Waveguide
        ax.plot( self.smoothedWG.x, self.smoothedWG.y, color="blue")
        unSmoothedWG = self.wpd.get("waveguide")
        ax.plot( unSmoothedWG.x, unSmoothedWG.y, color="blue", marker="x")
        ax.set_xlabel("$z$ (um)")
        ax.set_ylabel("$x$ (um)")
        fname = self.figFolder+"/"+self.overViewName
        fig.savefig(fname, bbox_inches="tight")
        print ("Figure written to %s"%(fname))

    def computeAngles( self ):
        if ( self.wpd is None ):
            print ("No WebPlotDigitizer object given!")
            return

        if ( self.smoothedWG is None ):
            self.createWG()

        self.angles = []
        self.angleX = []
        for dset in self.wpd.dsets:
            if ( dset.name == "waveguide" ):
                continue

            slope, interscept, pvalue, rvale, stderr = stats.linregress( dset.x, dset.y )
            closestIndx = np.argmin( np.abs(self.wgSlope.x - dset.x[0]) )
            self.angles.append( np.arctan(slope) - np.arctan(self.wgSlope.y[closestIndx]))
            self.angles[-1] *= 180.0/np.pi
            self.angleX.append( dset.x[0] )
        self.angles = np.array( self.angles )
        self.angleX = np.array( self.angleX )

    def plotWGSlope( self ):
        if ( self.wpd is None ):
            print ("No WebPlotDigitizer object given!")
            return

        if ( self.smoothedWG is None ):
            self.createWG()

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot( self.wgSlope.x, self.wgSlope.y, color="black")
        ax.set_xlabel("$z$ (um)")
        ax.set_ylabel("$dx/dz$")
        fname = self.figFolder + "/" + "wgSlope.pdf"
        fig.savefig(fname, bbox_inches="tight")
        print ("Figure written to %s"%(fname))

    def plotAngles( self ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot( self.angleX, self.angles, color="black")
        ax.set_xlabel("$z$ (um)")
        ax.set_ylabel("Ray angle (deg)")
        fname = self.figFolder + "/" + self.angleFname
        fig.savefig( fname, bbox_inches="tight")
        print ("Figure written to %s"%(fname))
