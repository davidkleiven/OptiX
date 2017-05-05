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
        self.limits = Limits()
        self.uid = 0
        self.name = ""
        self.xcrdLab = ""
        self.ycrdLab = ""
        self.quantityLab = ""
        self.cbLoc = None
        self.cbTick = None
        self.cbLog = True
        self.equalXYDim = True

    def getXlim( self ):
        return self.limits.xmin, self.limits.xmax

    def getYlim( self ):
        return self.limits.ymin, self.limits.ymax

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

    def projectionXY( self, cg ):
        maxval = np.max( self.data )
        minval = 1E-14*maxval
        if ( np.min(self.data) > minval ):
            minval = np.min(self.data)
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        xmin, xmax = self.getXlim()
        ymin, ymax = self.getYlim()

        extent = [xmin, xmax, ymin, ymax]
        aspect = self.aspectRatio( extent )
        if ( self.cbLog ):
            im = ax.imshow( self.data, extent=extent, aspect=aspect, cmap=self.cmap, norm=colors.LogNorm(minval,maxval) )
        else:
            im = ax.imshow( self.data, extent=extent, aspect=aspect, cmap=self.cmap )
        cb = fig.colorbar( im )

        ax.set_xlabel( self.xcrdLab )
        ax.set_ylabel( self.ycrdLab )
        #ax.set_xlabel( "$q_x \\backslash SI{}{\\backslash per\\backslash angstrom}$" )
        #ax.set_ylabel( "$q_y \\ SI{}{\\ per\\ angstrom}$" )
        cg.attach( fig, ax, "Figures/%sXY%d.svg"%(self.name, self.uid) )
        return cg

    def projectionX( self, cg, color="black", label="" ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        xmin, xmax = self.getXlim()
        x = np.linspace( xmin, xmax, len(self.data[0,:]) )
        ax.plot( x, self.data[int(self.center[0]),:], color=color, label=label)
        ax.spines["right"].set_visible( False )
        ax.spines["top"].set_visible( False )
        #ax.set_xlabel( "$q_x \\SI{}{\\per\\angstrom}$" )
        ax.set_xlabel( self.xcrdLab )
        ax.set_ylabel( self.quantityLab )
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        cg.attach( fig, ax, "Figures/%sX%d.svg"%( self.name, self.uid) )
        return cg

    def projectionY( self, cg, color="black", label="" ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ymin, ymax = self.getYlim()
        y = np.linspace( ymin, ymax, len(self.data[:,0]) )
        values = self.data[:,int(self.center[1])]
        values = values[::-1]
        ax.plot( y, values, color=color, label=label )
        ax.spines["right"].set_visible( False )
        ax.spines["top"].set_visible( False )
        #ax.set_xlabel( "$q_y \\SI{}{\\per\\angstrom}$" )
        ax.set_xlabel( self.ycrdLab )
        ax.set_ylabel( self.quantityLab )
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        cg.attach( fig, ax, "Figures/%sY%d.svg"%( self.name, self.uid) )
        return cg

class FarField(Plotter3D):
    def __init__( self ):
        Plotter3D.__init__( self )
        self.xcrdLab = r"\$q_x (\SI{}{\per\nano\meter})\$"
        self.ycrdlab = r"\$q_y (\SI{}{\per\nano\meter})\$"
        self.name = "FarField"
        self.quantityLab = "Intensity (a.u.)"

    def projectionXY( self, cg ):
        self.cbLoc = [np.min(self.data), np.max(self.data)]
        self.cbTick = ["$10^{%d}$"%(int(np.log10(np.min(self.data)))), "$10^{%d}$"%(int(np.log10(np.max(self.data))))]
        return Plotter3D.projectionXY( self, cg )

    def getXlim( self ):
        return self.limits.qmin, self.limits.qmax

    def getYlim( self ):
        return self.limits.qmin, self.limits.qmax

class Phase(Plotter3D):
    def __init__( self ):
        Plotter3D.__init__(self)
        #self.cbLoc = [-np.pi, 0.0, np.pi]
        #self.cbTick = ["-$\pi$", "$0$", "$\pi$"]
        self.xcrdLab =  "$x \\backslash SI{}{\\backslash micro\\backslash meter}$"
        self.ycrdlab = "$y \\backslash SI{}{\\backslash micro\\backslash meter}$"
        self.name = "phase"
        self.cbLog = False
        self.quantityLab = "Phase (rad)"

class IntensityPlot(Plotter3D):
    def __init_( self ):
        Plotter3D.__init__( self )
        self.xcrdLab =  r"\$x \SI{}{\micro\meter}\$"
        self.ycrdlab = "\$y \SI{}{\micro\meter}\$"
        self.name = "intensity"
        self.quantityLab = "Intensity (a.u.)"
        self.cbLog = False
        print (self.cbLog)
