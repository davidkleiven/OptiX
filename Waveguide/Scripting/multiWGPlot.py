import numpy as np
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
            minval = 1E-8*maxval
            im = ax.imshow( intensity, extent=extent, origin="lower", cmap=self.cmap, norm=mpl.colors.LogNorm(minval,maxval) )
        else:
            im = ax.imshow( intensity, extent=extent, origin="lower", cmap=self.cmap )

        ax.set_aspect( (extent[1]-extent[0])/(extent[3]-extent[2]) )
        ax.autoscale(False)
        ax.set_xlabel("Arclength (\$\SI{}{\milli\meter}\$)")
        ax.set_ylabel("Transverse position (\$\SI{}{\\nano\meter}\$)")
        fig.colorbar( im )
        fig, ax = self.addMiniatyrGeometry( fig, ax )
        return fig, self.drawSeparationLines( ax, "white" )

    def plotTransmittivity(self, trans ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        z = np.linspace(self.zmin*self.mm, self.zmax*self.mm, len(trans))
        ax.plot( z, np.log(trans), color="black", lw=2)
        ax.set_xlabel("Arclength (\$\SI{}{\milli\meter}\$)")
        ax.set_ylabel("ln( Transmittivity )")
        return fig, self.drawSeparationLines( ax, "#2b8cbe")

    def addMiniatyrGeometry( self, fig, ax ):
        if ( not path.exists( self.miniatyrGeo ) ):
            print ("Could not find image file containing the geometry!")
            return
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        im = plt.imread( self.miniatyrGeo, format="png" )
        #ax.imshow( im )
        imgbox = offsetbox.OffsetImage( im, zoom=0.002)
        newax = inl.inset_axes( ax, width="40%", height="40%", loc=1)
        newax.imshow( im )
        newax.axis("off")
        #ab = offsetbox.AnnotationBbox( imgbox, )
        #ax.add_artist( imgbox )
        #width = 0.3*(xmax-xmin)
        #height = 0.25*(ymax-ymin)
        #rect = [xmax-width, ymax.height, width, height]


        #newax = fig.axes( [xmax-width, ymax-height, width, height] )
        #newax.imshow( im )
        #newax.axes("off")
        return fig, ax
