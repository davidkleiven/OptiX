import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import offsetbox
from os import path
from mpl_toolkits.axes_grid import inset_locator as inl

class MultiWGPlot:
    def __init__( self ):
        self.R = []
        self.angles = []
        self.xmin = 0.0
        self.xmax = 1.0
        self.zmin = 0.0
        self.zmax = 1.0
        self.cmap = "nipy_spectral"
        self.miniatyrGeo = ""
        self.miniatyr = True

        # Define some units
        self.nm = 1.0
        self.um = 1E-3
        self.mm = 1E-6

    def drawSeparationLines(self, ax, color ):
        # Set separation lines
        assert( len(self.R) == len(self.angles))
        endpos = 0.0
        for i in range(0, len(self.R)-1):
            endpos += self.R[i]*self.angles[i]*np.pi/180.0
            ax.axvline( endpos, ls="--", color=color, lw=2)
        return ax

    def plotIntensity(self, intensity, scale="linear" ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        extent = [self.zmin*self.mm, self.zmax*self.mm, self.xmin, self.xmax]
        if ( scale == "log" ):
            maxval = np.max(intensity)
            minval = 1E-5*maxval
            im = ax.imshow( intensity, extent=extent, aspect="auto", origin="lower", cmap=self.cmap, norm=mpl.colors.LogNorm(minval,maxval) )
        else:
            im = ax.imshow( intensity, extent=extent, aspect="auto", origin="lower", cmap=self.cmap )

        ax.set_aspect( (extent[1]-extent[0])/(extent[3]-extent[2]) )
        ax.autoscale(False)
        ax.set_xlabel("Arclength (\$\SI{}{\milli\meter}\$)")
        ax.set_ylabel("Transverse position (\$\SI{}{\\nano\meter}\$)")
        fig.colorbar( im )
        if ( self.miniatyr ):
            fig = self.addMiniatyrGeometry( fig, ax )
        return fig, self.drawSeparationLines( ax, "white" )

    def plotTransmittivity(self, trans ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        z = np.linspace(self.zmin*self.mm, self.zmax*self.mm, len(trans))
        ax.plot( z, np.log(trans), color="black", lw=2)
        ax.set_xlabel("Arclength (\$\SI{}{\milli\meter}\$)")
        ax.set_ylabel("ln( Transmittivity )")
        ax.spines["right"].set_visible( False )
        ax.spines["top"].set_visible( False )
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        if ( self.miniatyr ):
            fig = self.addMiniatyrGeometry( fig, ax )
        return fig, self.drawSeparationLines( ax, "#2b8cbe")

    def addMiniatyrGeometry( self, fig, ax ):
        if ( not path.exists( self.miniatyrGeo ) ):
            print ("Could not find image file containing the geometry!")
            return
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        im = plt.imread( self.miniatyrGeo, format="png" )
        width = 0.5
        height = 0.4
        rect = [0.4, 0.55, 0.35, 0.35]
        newax = fig.add_axes( rect, frameon=False )
        newax.imshow( im )
        newax.axis("off")
        return fig
