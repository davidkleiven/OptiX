import sys
from scipy import ndimage as nd
from matplotlib import pyplot as plt

class Limits:
    def __init__( self ):
        self.xmin = 0.0
        self.ymin = 0.0
        self.zmin = 0.0

        self.xmax = 1.0
        self.ymax = 1.0
        self.zmax = 1.0

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
