import sys
from scipy import ndimage as nd
from matplotlib import pyplot as plt
from matplotlib import colors
import numpy as np

class Limits:
    def __init__( self ):
        self.xmin = 0.0
        self.ymin = 0.0
        self.zmin = 0.0

        self.xmax = 1.0
        self.ymax = 1.0
        self.zmax = 1.0

        self.qmin = 0.0
        self.qmax = 0.0

class Plotter3D:
    def __init__( self ):
        self.data = None
        self.center = [0,0,0]
        self.cmap = "nipy_spectral"
        self.limits = None
        self.uid = 0
        self.name = ""

    def centerOfMass( self ):
        self.center = nd.measurements.center_of_mass( self.data )

    def aspectRatio( self, extent ):
        return ( extent[1] - extent[0] )/( extent[3] - extent[2])

    def sliceZ( self, cg ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        extent = [self.limits.xmin/1000.0, self.limits.xmax/1000.0, self.limits.ymin/1000.0, self.limits.ymax/1000.0]
        aspect = self.aspectRatio( extent )
        ax.imshow( self.data[:,:, self.center[2]], extent=extent, aspect=aspect, cmap=self.cmap )
        ax.set_xlabel( "$x \SI{}{\micro\meter}$" )
        ax.set_ylabel( "$y \SI{}{\micro\meter}$" )
        cg.attach( fig, ax, "Figures/%ssliceZ%d.svg"%(self.name, self.uid) )
        return cg

    def farField( self, cg ):
        maxval = np.max( self.data )
        minval = 1E-8*maxval
        if ( np.min(self.data) > minval ):
            minval = np.min(self.data)
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        extent = [self.limits.qmin*10.0, self.limits.qmax*10.0, self.limits.qmin*10.0, self.limits.qmax*10.0]
        aspect = self.aspectRatio( extent )
        im = ax.imshow( self.data, extent=extent, aspect=aspect, cmap=self.cmap, norm=colors.LogNorm(minval,maxval) )
        cb = fig.colorbar( im )
        labels = ["$10^{%d}$"%(int(np.log10(minval))), "$10^{%d}$"%(int(np.log10(maxval)))]
        loc = [minval, maxval]
        cb.set_ticks(loc)
        cb.set_ticklabels(labels)
        ax.set_xlabel( "$q_x \backslash SI{}{\backslash per\backslash angstrom}$" )
        #ax.set_ylabel( "$q_y \\ SI{}{\\ per\\ angstrom}$" )
        cg.attach( fig, ax, "Figures/%s%d.svg"%(self.name, self.uid) )
        return cg

    def projectionX( self, cg ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        q = np.linspace( self.limits.qmin*10.0, self.limits.qmax*10.0, len(self.data[0,:]) )
        ax.plot( q, self.data[int(self.center[0]),:], color="black")
        ax.spines["right"].set_visible( False )
        ax.spines["top"].set_visible( False )
        #ax.set_xlabel( "$q_x \\SI{}{\\per\\angstrom}$" )
        cg.attach( fig, ax, "Figures/%sQx%d.svg"%( self.name, self.uid) )
        return cg

    def projectionY( self, cg ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        q = np.linspace( self.limits.qmin*10.0, self.limits.qmax*10.0, len(self.data[:,0]) )
        ax.plot( q, self.data[:,int(self.center[1])], color="black")
        ax.spines["right"].set_visible( False )
        ax.spines["top"].set_visible( False )
        #ax.set_xlabel( "$q_y \\SI{}{\\per\\angstrom}$" )
        cg.attach( fig, ax, "Figures/%sQy%d.svg"%( self.name, self.uid) )
        return cg
