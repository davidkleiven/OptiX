import sys
import numpy as np
from matplotlib import patches as ptc
from matplotlib import pyplot as plt
import h5py as h5

class Geometry:
    def __init__( self, hdf ):
        if ( "geometrySourceVolume" in hdf.keys() ):
            self.srcPlane = np.array( hdf.get("geometrySourceVolume") )
        else:
            print ("Source plane not in HDF5 file!")
            self.srcPlane = None

        if ( "boxGeo" in hdf.keys() ):
            self.boxGeo = np.array( hdf.get("boxGeo") )
        else:
            print ("Box geometry not in HDF5 file!")
            self.boxGeo = None

        if ( "refPGeo" in hdf.keys() ):
            self.refPlane = np.array( hdf.get("refPGeo") )
        else:
            print ("Reflection plane not in HDf5 file!")
            self.refPlane = None

        if ( "trGeo" in hdf.keys() ):
            self.trPlane = np.array( hdf.get("trPGeo") )
        else:
            print ("Transmission plane not in HDF5 file!")
            self.trPlane = None

        if ( "eps" in hdf.keys() ):
            self.eps = np.array( hdf.get("eps") )
        else:
            print ("Epsilon not in HDF5 file!")
            self.eps = None

        # Initialize figure
        self.fig = plt.figure()
        self.ax1 = self.fig.add_subplot(2,2,1)
        self.ax2 = self.fig.add_subplot(2,2,2)
        self.ax3 = self.fig.add_subplot(2,2,3)

    def plotEps( self ):
        if ( self.eps is None ):
            return
        assert( len(self.eps.shape) == 3 )
        projXY = np.sum( self.eps , axis=2 )
        self.ax1.imshow( projXY, cmap="bone")#, aspect="auto" )

        projXZ = np.sum( self.eps, axis=1 )
        self.ax2.imshow( projXZ, cmap="bone")#, aspect="auto")

        projYZ = np.sum( self.eps, axis=0 )
        self.ax3.imshow( projYZ, cmap="bone")#, aspect="auto")

    def plotSrc( self ):
        if ( self.srcPlane is None ):
            return
        x = [self.srcPlane[0], self.srcPlane[3]]
        y = [self.srcPlane[1], self.srcPlane[4]]
        z = [self.srcPlane[2], self.srcPlane[5]]
        color="#e41a1c"
        self.ax1.plot( y, x, color=color, lw=4, label="Source" )
        self.ax2.plot( z, x, color=color, lw=4, label="Source" )
        self.ax3.plot( z, y, color=color, lw=4, label="Source" )

    def plotRefPlane( self ):
        if ( self.refPlane is None ):
            return
        x = [self.refPlane[0], self.refPlane[3]]
        y = [self.refPlane[1], self.refPlane[4]]
        z = [self.refPlane[2], self.refPlane[5]]
        color="#377eb8"
        self.ax1.plot( y, x, color=color, lw=3, label="RefFlux" )
        self.ax2.plot( z, x, color=color, lw=3, label="RefFlux" )
        self.ax3.plot( z, y, color=color, lw=3, label="RefFlux" )

    def plotTransPlane( self ):
        if ( self.trPlane is None ):
            return
        x = [self.trPlane[0], self.trPlane[3]]
        y = [self.trPlane[1], self.trPlane[4]]
        z = [self.trPlane[2], self.trPlane[5]]
        color="#4daf4a"
        self.ax1.plot( y, x, color=color, lw=2, label="TransFlux" )
        self.ax2.plot( z, x, color=color, lw=2, label="TransFlux" )
        self.ax3.plot( z, y, color=color, lw=2, label="TransFlux" )

    def plot( self ):
        self.plotEps()
        self.plotSrc()
        self.plotRefPlane()
        self.plotTransPlane()
        self.ax3.legend(fancybox=True, shadow=True, bbox_to_anchor=(3.0,1.0))

    def plotFluxBox( self ):
        if ( self.boxGeo is None ):
            return
        color = "#984ea3"
        rect = ptc.Rectangle( (self.boxGeo[0],self.boxGeo[1]), self.boxGeo[3]-self.boxGeo[0], self.boxGeo[4]-self.boxGeo[1], fill=False, ec=color )
        self.ax1.add_patch( rect )

        rect = ptc.Rectangle( (self.boxGeo[2],self.boxGeo[0]), self.boxGeo[5]-self.boxGeo[2], self.boxGeo[3]-self.boxGeo[0], fill=False, ec=color )
        self.ax2.add_patch( rect )

        rect = ptc.Rectangle( (self.boxGeo[2],self.boxGeo[1]), self.boxGeo[5]-self.boxGeo[2], self.boxGeo[4]-self.boxGeo[1], fill=False, ec=color )
        self.ax2.add_patch( rect )

def main( argv ):
    if ( len(argv) != 1 ):
        print ("Usage: python geometry.py <hdf5 file>" )
        return

    with h5.File(argv[0], 'r') as hf:
        geo = Geometry( hf )

    geo.plot()
    plt.show()

if __name__ == "__main__":
    main( sys.argv[1:] )
