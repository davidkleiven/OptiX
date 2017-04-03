import sys
import numpy as np
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 18
mpl.rcParams["axes.unicode_minus"] = False
from matplotlib import pyplot as plt
import h5py as h5
import subprocess

class Polarization:
    def __init__( self ):
        self.I = []
        self.Q = []
        self.U = []
        self.V = []
        self.wavelength = []
        self.theta = []

    def load( self, hfile ):
        key = "StokesI0"
        indx = 0
        self.theta = np.array( hfile.get("SrTheta"))
        self.I = np.array( hfile.get("stokesIAvg") )
        self.Q = np.array( hfile.get("stokesQAvg") )
        self.U = np.array( hfile.get("stokesUAvg"))
        self.V = np.array( hfile.get("stokesVAvg") )
        freqs = np.array( hfile.get("spectrumFreqs"))
        freq = np.linspace(freqs[0], freqs[0]+freqs[1]*freqs[2], self.I.shape[1])
        vxsize = np.array( hfile.get("vxSizeNM"))[0]
        self.wavelength = vxsize/freq

    def plot( self, freq ):
        if ( indx >= len(self.I) ):
            raise Exception("Index out of bounds in plot function...")
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot( self.theta, self.Q[:,freq]/self.I[:,freq], label="Q")
        ax.plot( self.theta, self.U[:,freq]/self.I[:,freq], label="U")
        ax.plot( self.theta, self.V[:,freq]/self.I[:,freq], label="V")
        ax.legend( loc="best", frameon=False, labelspacing=0.05)
        return fig, ax

    def plotDistributions( self ):
        fig = plt.figure()
        ax = []
        ax.append( fig.add_subplot(1,3,1) )
        WAVEL, THETA = np.meshgrid( self.wavelength, self.theta )
        THETA *= 180.0/np.pi
        maxval = np.max( [(self.Q/self.I).max(),(self.U/self.I).max(),(self.V/self.I).max()] )
        minval = np.min( [(self.Q/self.I).min(), (self.U/self.I).min(), (self.V/self.I).min()] )
        cmap="RdYlGn"
        cmap="Spectral_r"
        ncol = 255
        ax[0].contourf( THETA, WAVEL, self.Q/self.I, ncol, cmap=cmap, vmin=minval,vmax=maxval)
        ax[0].set_ylabel("Wavelength (nm)")
        ax[0].set_xlabel("Scattering angle (deg)")
        ax[0].set_title("Q")

        ax.append( fig.add_subplot(1,3,2) )
        ax[1].contourf( THETA, WAVEL, self.U/self.I, ncol, cmap=cmap, vmin=minval,vmax=maxval)
        ax[1].set_yticks([])
        ax[1].set_title("U")

        ax.append( fig.add_subplot(1,3,3) )
        im = ax[2].contourf( THETA, WAVEL, self.V/self.I, ncol, cmap=cmap, vmin=minval,vmax=maxval)
        ax[2].set_yticks([])
        ax[2].set_title("V")
        fig.colorbar(im)

    def freqAvg( self ):
        return self.I.sum(axis=1), self.Q.sum(axis=1), self.U.sum(axis=1), self.V.sum(axis=1)

    def plotFrequencyAverage( self ):
        I,Q,U,V = self.freqAvg()
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot( self.theta*180.0/np.pi, Q/I, color="#e41a1c", label="Q")
        ax.plot( self.theta*180.0/np.pi, U/I, color="#377eb8", label="U")
        ax.plot( self.theta*180.0/np.pi, V/I, color="#4daf4a", label="V")
        ax.set_xlabel("Scattering angle (deg)")
        ax.set_ylabel("Stokes parameter (Q,U,V)")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.xaxis.set_ticks_position("bottom")
        ax.yaxis.set_ticks_position("left")
        ax.legend(loc="best", frameon=False, labelspacing=0.05)

    def polarizationEllipse( self, ang, freq=None ):
        if ( freq is None ):
            I,Q,U,V = self.freqAvg()
        else:
            I = self.I[:,freq]
            Q = self.Q[:,freq]
            U = self.U[:,freq]
            V = self.V[:,freq]
        psi = 0.5*np.arctan2( U[ang], Q[ang] )
        xsi = 0.5*np.arctan2( V[ang], np.sqrt(U[ang]**2 + Q[ang]**2) )
        center = (0.0,0.0)
        height = np.tan(xsi)
        width = 2.0
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        maxval = np.max([height,width])/2.0
        ax.set_xlim(-maxval,maxval)
        ax.set_ylim(-maxval,maxval)
        ellipse = mpl.patches.Ellipse( center, width, height, angle=psi*180.0/np.pi, edgecolor="black", fc="none", lw=2 )
        ax.add_patch(ellipse)
        ax.set_title("Scattering angle: %.2f"%(self.theta[ang]*180.0/np.pi))
        return fig, ax

    def animationPolarizationEllipse( self, freq=None ):
        dirname = "tempAnimDir"
        subprocess.call(["mkdir", "-p", dirname])
        for i in range(0, len(self.theta)):
            fig, ax = self.polarizationEllipse(i, freq=freq)
            fig.savefig(dirname+"/img%d.png"%(i))
            plt.close()

        # Create animation
        subprocess.call(["ffmpeg", "-r", "2", "-i", dirname+"/img%d.png", "-r", "10", "data/polarizationEllipse.ogg"])
        subprocess.call(["rm", "-r", dirname])


def main( argv ):
    if ( len(argv) < 1 ) or (len(argv)>2 ):
        print ("Usage: python polarization.py hfilfe.h5")
        return 1

    createAnimation = False
    for arg in argv:
        if ( arg.find("--anim") != -1 ):
            createAnimation = True

    polPlots = Polarization()
    with h5.File( argv[0], 'r' ) as hf:
        polPlots.load(hf)

    try:
        polPlots.plotDistributions()
        polPlots.plotFrequencyAverage()
        polPlots.polarizationEllipse( 0 )
        if ( createAnimation ):
            polPlots.animationPolarizationEllipse(freq=75)
        plt.show()
    except Exception as exc:
        print (str(exc))
        return 1
    return 0

if __name__ == "__main__":
    main( sys.argv[1:] )
