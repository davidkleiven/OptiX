import numpy as np
from matplotlib import pyplot as plt

class MultiWGPlot:
    def __init__( self ):
        self.R = []
        self.angles = []
        self.xmin = 0.0
        self.xmax = 1.0
        self.zmin = 0.0
        self.zmax = 1.0
        self.cmap = "nipy_spectral"

        # Define some units
        self.nm = 1.0
        self.um = 1E-3
        self.mm = 1E-6

    def drawSeparationLines(self, ax, color ):
        # Set separation lines
        assert( len(self.R) == len(self.angles))
        for i in range(0, len(self.R)-1):
            endpos = self.R[i]*self.mm*self.angles[i]*np.pi/180.0
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
        return fig, self.drawSeparationLines( ax, "white" )

    def plotTransmittivity(self, trans ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        z = np.linspace(self.zmin*self.mm, self.zmax*self.mm, len(trans))
        ax.plot( z, np.log(trans), color="black", lw=2)
        ax.set_xlabel("Arclength (\$\SI{}{\milli\meter}\$)")
        ax.set_ylabel("ln( Transmittivity )")
        return fig, self.drawSeparationLines( ax, "#2b8cbe")
