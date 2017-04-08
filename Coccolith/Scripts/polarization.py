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

        self.IphiTheta = [] # I as a function of phi and theta
        self.QphiTheta = [] # Q as a function of phi and theta
        self.UphiTheta = [] # U as a function of phi and theta
        self.VphiTheta = [] # V as a function of phi and theta
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
        ax.plot( self.theta, self.I[:,freq], label="I")
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

    def readPhiAndTheta( self, hfile, freqIndx ):
        # Count how many phi angles there are
        nphi = 0
        key = "stokesI%d"%(nphi)
        ntheta = 0
        while ( key in hfile.keys() ):
            if ( nphi == 0 ):
                ntheta = np.array( hfile.get(key) ).shape[0]
            nphi += 1
            key = "stokesI%d"%(nphi)

        # Initialize arrays
        self.IphiTheta = np.zeros((ntheta,nphi))
        self.QphiTheta = np.zeros((ntheta,nphi))
        self.UphiTheta = np.zeros((ntheta,nphi))
        self.VphiTheta = np.zeros((ntheta,nphi))

        # Read arrays
        for i in range(0, nphi):
            self.IphiTheta[:,i] = np.array( hfile.get("stokesI%d"%(i)) )[:,freqIndx]
            self.QphiTheta[:,i] = np.array( hfile.get("stokesQ%d"%(i)) )[:,freqIndx]
            self.UphiTheta[:,i] = np.array( hfile.get("stokesU%d"%(i)) )[:,freqIndx]
            self.VphiTheta[:,i] = np.array( hfile.get("stokesV%d"%(i)) )[:,freqIndx]

    def plotPolarizationPhiTheta( self ):
        fig = plt.figure()
        axI = fig.add_subplot(2,2,1)
        axQ = fig.add_subplot(2,2,2)
        axU = fig.add_subplot(2,2,3)
        axV = fig.add_subplot(2,2,4)
        phi = np.linspace(0,2.0*np.pi,self.IphiTheta.shape[1])
        PHI, THETA = np.meshgrid(phi, self.theta)
        cmap ="RdYlBu_r"
        cmap = "Spectral"
        figdebug = plt.figure()
        axdebug = figdebug.add_subplot(1,1,1)
        inorm = self.IphiTheta/self.IphiTheta.max()
        axdebug.contourf( THETA, PHI, inorm, 512, cmap="inferno")#, norm=mpl.colors.LogNorm())
        extent = [0.0,180,0.0,360]
        imI = axI.imshow( inorm.T[:,::-1], cmap="inferno", extent=extent, aspect="auto")
        axI.set_title("I")
        #fig.colorbar(imI)
        imQ = axQ.imshow( self.QphiTheta.T[:,::-1]/self.IphiTheta.T[:,::-1], extent=extent, cmap=cmap, aspect="auto" )
        #imQ = axQ.imshow( self.QphiTheta.T[:,::-1], extent=extent, cmap=cmap, aspect="auto" )
        axQ.set_title("Q")
        #fig.colorbar(imQ)
        imU = axU.imshow( self.UphiTheta.T[:,::-1]/self.IphiTheta.T[:,::-1], extent=extent, cmap=cmap, aspect="auto" )
        axU.set_title("U")
        #fig.colorbar(imQ)
        imV = axV.imshow( self.VphiTheta.T[:,::-1]/self.IphiTheta.T[:,::-1], extent=extent, cmap=cmap, aspect="auto")
        axV.set_title("V")
        #fig.colorbar(imV)
        axV.set_xlabel("Scattering angle, \$\\theta\$ (deg)")
        axU.set_ylabel("Azimuthal angle, \$\phi\$ (deg)")

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
        polPlots.readPhiAndTheta( hf, 2 )

    try:
        polPlots.plotDistributions()
        polPlots.plotFrequencyAverage()
        polPlots.polarizationEllipse( 0 )
        polPlots.plotPolarizationPhiTheta()
        if ( createAnimation ):
            polPlots.animationPolarizationEllipse(freq=1)
        plt.show()
    except Exception as exc:
        print (str(exc))
        return 1
    return 0

if __name__ == "__main__":
    main( sys.argv[1:] )
