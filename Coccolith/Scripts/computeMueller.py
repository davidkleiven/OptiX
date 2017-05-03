import sys
sys.path.append("Scripts")
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import polarization as plz
import h5py as h5
from scipy import interpolate as interp

class IncidentStokes:
    def __init__( self ):
        self.IInc = []
        self.QInc = []
        self.UInc = []
        self.VInc = []

class MuellerMatrix:
    def __init__( self ):
        self.pol = []
        self.stokesInc = []
        self.theta = []
        self.phi = []

    def load( self, hfile ):
        counter = 0
        self.theta = np.array( hfile.get("theta") )
        self.phi = np.array( hfile.get("phi") )
        key = "Run%d"%(counter)
        while (key in hfile.keys() ):
            pol = plz.Polarization()
            pol.readPhiAndTheta( hfile.get(key), 1 )
            self.pol.append(pol)
            incStk = IncidentStokes()
            incStk.Iinc = np.array( hfile.get(key).get("StokesIInc") )
            incStk.Qinc = np.array( hfile.get(key).get("StokesQInc") )
            incStk.Uinc = np.array( hfile.get(key).get("StokesUInc") )
            incStk.Vinc = np.array( hfile.get(key).get("StokesVInc") )
            self.stokesInc.append(incStk)
            counter += 1
            key = "Run%d"%(counter)

    def getIncStokesMatrix( self, phiIndx ):
        assert( len(self.stokesInc) == 4 )
        mat = np.zeros((4,4))
        for i in range(0,4):
            mat[0,i] = self.stokesInc[i].Iinc[phiIndx]
            mat[1,i] = self.stokesInc[i].Qinc[phiIndx]
            mat[2,i] = self.stokesInc[i].Uinc[phiIndx]
            mat[3,i] = self.stokesInc[i].Vinc[phiIndx]
        return mat

    def getScatteredStokesMatrix( self, thetaIndx, phiIndx ):
        assert( len(self.pol) == 4 )
        mat = np.zeros((4,4))
        for i in range(0,4):
            I = self.pol[i].IphiTheta[thetaIndx,phiIndx]
            I = 1.0
            mat[0,i] = self.pol[i].IphiTheta[thetaIndx,phiIndx]/I
            mat[1,i] = self.pol[i].QphiTheta[thetaIndx,phiIndx]/I
            mat[2,i] = self.pol[i].UphiTheta[thetaIndx,phiIndx]/I
            mat[3,i] = self.pol[i].VphiTheta[thetaIndx,phiIndx]/I
        return mat

    def mueller( self, thetaIndx, phiIndx ):
        Sinc = self.getIncStokesMatrix(phiIndx)
        Sscat = self.getScatteredStokesMatrix( thetaIndx, phiIndx )
        return Sscat.dot( np.linalg.inv(Sinc) )

    def projectionMatrix( self, thetaIndx, phiIndx ):
        Sinc = self.getIncStokesMatrix(phiIndx)
        Sscat = self.getScatteredStokesMatrix( thetaIndx, phiIndx )
        for i in range(0,4):
            Sscat[:,i] /= Sscat[0,i]
        return 0.5*(Sscat.T).dot( Sinc )

    def plot( self ):
        ntheta = self.pol[0].IphiTheta.shape[0]
        nphi = self.pol[0].IphiTheta.shape[1]

        muellerAng = [[],[],[],[]]
        for i in range(0,4):
            for j in range(0,4):
                muellerAng[i].append(np.zeros((ntheta,nphi)))

        for theta in range(0,ntheta):
            for phi in range(0,nphi):
                mat = self.mueller(theta,phi)
                mat /= mat[0,0]
                for i in range(0,4):
                    for j in range(0,4):
                        muellerAng[i][j][theta,phi] = mat[i,j]

        fig = plt.figure()
        axes = []
        for i in range(0,16):
            axes.append( fig.add_subplot(4,4,i+1) )

        cmap = "Spectral_r"
        #cmap = "coolwarm"
        maxval = -1E30
        minval = 1E30
        # Locate max and min value for colorbar scaling
        for i in range(0,4):
            for j in range(0,4):
                if ( muellerAng[i][j].max() > maxval ):
                    maxval = muellerAng[i][j].max()
                if ( muellerAng[i][j].min() < minval ):
                    minval = muellerAng[i][j].min()

        THETA, PHI = np.meshgrid( self.theta, self.phi )
        for i in range(0,4):
            for j in range(0,4):
                rectInterp = interp.RectBivariateSpline( self.theta[::-1], self.phi, muellerAng[i][j][::-1,:] )
                imgsize = 128
                thetaInterp = np.linspace( self.theta.min(), self.theta.max(), imgsize )
                phiInterp = np.linspace( self.phi.min(), self.phi.max(), imgsize )
                img = np.zeros((imgsize,imgsize))
                for th in range(0,len(thetaInterp)):
                    for ph in range(0,len(phiInterp)):
                        img[ph,th] = rectInterp(thetaInterp[th],phiInterp[ph])
                im = axes[i*4+j].imshow( img, cmap=cmap, aspect="auto", vmin=minval, vmax=maxval)
                #im = axes[i*4+j].contourf( THETA, PHI, muellerAng[i][j].T, 256, cmap=cmap, vmin=minval, vmax=maxval)
        for i in range(0,16):
            axes[i].xaxis.set_ticks([])
            axes[i].yaxis.set_ticks([])

        cbar_ax = fig.add_axes([0.15, 0.95, 0.7, 0.05])
        fig.colorbar( im, orientation="horizontal", cax=cbar_ax )

    def plotStokesVectorVSPhi( self, runIndx, thetaIndx ):
        fig = plt.figure()
        ax1 = fig.add_subplot(1,2,1)
        ax1.plot( self.phi, self.pol[runIndx].IphiTheta[thetaIndx,:], label="I")
        ax1.plot( self.phi, self.pol[runIndx].QphiTheta[thetaIndx,:], label="Q")
        ax1.plot( self.phi, self.pol[runIndx].UphiTheta[thetaIndx,:], label="U")
        ax1.plot( self.phi, self.pol[runIndx].VphiTheta[thetaIndx,:], label="V")
        ax1.legend()

        ax2 = fig.add_subplot(1,2,2)
        ax2.plot( self.stokesInc[runIndx].Iinc, label="I")
        ax2.plot( self.stokesInc[runIndx].Qinc, label="Q")
        ax2.plot( self.stokesInc[runIndx].Uinc, label="U")
        ax2.plot( self.stokesInc[runIndx].Vinc, label="V")
        ax2.legend()


    def degreeOfLinearPolarizationUnpolarized( self, phi ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]
        for j in range(0,len(phi)):
            phiIndx = np.argmin( np.abs(self.phi-phi[j]) )
            degreeLinPol = np.zeros( len(self.theta) )
            for i in range(0,len(self.theta) ):
                mueller = self.mueller(i,phiIndx)
                degreeLinPol[i] = np.sqrt( mueller[1,0]**2 + mueller[2,0]**2 )
                degreeLinPol[i] /= mueller[0,0]

            ax.plot( self.theta*180.0/np.pi, degreeLinPol, color=colors[j%len(colors)],
            label="\$\SI{%d}{\degree}\$"%(phi[j]*180.0/np.pi) )
        ax.set_xlabel("Scattering angle, \$\\theta \$ (deg)")
        ax.set_ylabel( "\$\\frac{\sqrt{M_{21}^2+M_{31}^2}}{M_{11}}\$" )
        ax.legend( loc="best", frameon=False, labelspacing=0.05 )
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.xaxis.set_ticks_position("bottom")
        ax.yaxis.set_ticks_position("left")

    def degreeOfLinPolPhiAveraged( self ):
        degreeLinPol = np.zeros( len(self.theta) )
        for j in range(0,len(self.phi)):
            for i in range(0,len(self.theta) ):
                mueller = self.mueller(i,j)
                degreeLinPol[i] += np.sqrt( mueller[1,0]**2 + mueller[2,0]**2 )/mueller[0,0]
        degreeLinPol /= len(self.phi)
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot( self.theta*180.0/np.pi, degreeLinPol, color="black" )
        ax.set_xlabel("Scattering angle, \$\\theta \$ (deg)")
        ax.set_ylabel( "\$\Big\langle \\frac{\sqrt{M_{21}^2+M_{31}^2}}{M_{11}} \Big\\rangle_\phi\$" )
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.xaxis.set_ticks_position("bottom")
        ax.yaxis.set_ticks_position("left")

    def muellerMatrixSingleElements( self, phi ):
        fig = plt.figure()
        ax1 = fig.add_subplot(2,2,1)
        ax2 = fig.add_subplot(2,2,2)
        ax3 = fig.add_subplot(2,2,3)
        ax4 = fig.add_subplot(2,2,4)
        phiIndx = np.argmin( np.abs(self.phi-phi) )
        M11 = np.zeros(len(self.theta))
        M12 = np.zeros(len(self.theta))
        M13 = np.zeros(len(self.theta))
        M14 = np.zeros(len(self.theta))
        M21 = np.zeros(len(self.theta))
        M22 = np.zeros(len(self.theta))
        M23 = np.zeros(len(self.theta))
        M24 = np.zeros(len(self.theta))
        M31 = np.zeros(len(self.theta))
        M32 = np.zeros(len(self.theta))
        M33 = np.zeros(len(self.theta))
        M34 = np.zeros(len(self.theta))
        M41 = np.zeros(len(self.theta))
        M42 = np.zeros(len(self.theta))
        M43 = np.zeros(len(self.theta))
        M44 = np.zeros(len(self.theta))
        for i in range(0,len(self.theta) ):
            M = self.mueller(i,phiIndx)
            M11[i] = M[0,0]
            M12[i] = M[0,1]
            M13[i] = M[0,2]
            M14[i] = M[0,3]
            M21[i] = M[1,0]
            M22[i] = M[1,1]
            M23[i] = M[1,2]
            M24[i] = M[1,3]
            M31[i] = M[2,0]
            M32[i] = M[2,1]
            M33[i] = M[2,2]
            M34[i] = M[2,3]
            M41[i] = M[3,0]
            M42[i] = M[3,1]
            M43[i] = M[3,2]
            M44[i] = M[3,3]
        angle = self.theta*180.0/np.pi

        # Determine scale
        maxVals = np.max( [np.max(M11),np.max(M22),np.max(M33),np.max(M44)] )
        log10 = int( np.log10(maxVals) )
        if ( log10 < 0 ):
            log10 -= 1
        scale = 10**log10
        ax1.plot( angle, M11/scale, label="\$M_{11}\$", color="#e41a1c")
        ax1.plot( angle, M22/scale, label="\$M_{22}\$", color="#377eb8")
        ax1.plot( angle, M33/scale, label="\$M_{33}\$", color="#4daf4a")
        ax1.plot( angle, M44/scale, label="\$M_{44}\$", color="#ff7f00")
        ax1.set_xlabel("Scattering angle, \$\\theta\$ (deg)")
        ax1.set_ylabel("\$\\times 10^{%d}\$"%(log10) )
        ax1.legend(loc="best", frameon=False, labelspacing=0.05)
        ax1.spines["right"].set_visible(False)
        ax1.spines["top"].set_visible(False)
        ax1.xaxis.set_ticks_position("bottom")
        ax1.yaxis.set_ticks_position("left")

        maxVals = np.max( [np.max(M12),np.max(M21),np.max(M34),np.max(M43)] )
        log10 = int( np.log10(maxVals) )
        if ( log10 < 0 ):
            log10 -= 1
        scale = 10**log10
        ax2.plot( angle, M12/scale, label="\$M_{12}\$", color="#e41a1c" )
        ax2.plot( angle, M21/scale, label="\$M_{21}\$", color="#377eb8" )
        ax2.plot( angle, M34/scale, label="\$M_{34}\$", color="#4daf4a" )
        ax2.plot( angle, M43/scale, label="\$M_{43}\$", color="#ff7f00" )
        ax2.set_ylabel("\$\\times 10^{%d}\$"%(log10) )
        ax2.legend( loc="best", frameon=False, labelspacing=0.05 )
        ax2.spines["right"].set_visible(False)
        ax2.spines["top"].set_visible(False)
        ax2.xaxis.set_ticks_position("bottom")
        ax2.yaxis.set_ticks_position("left")

        maxVals = np.max( [np.max(M13),np.max(M31),np.max(M41),np.max(M14)] )
        log10 = int( np.log10(maxVals) )
        if ( log10 < 0 ):
            log10 -= 1
        scale = 10**log10
        ax3.plot( angle, M13/scale, label="\$M_{13}\$", color="#e41a1c" )
        ax3.plot( angle, M31/scale, label="\$M_{31}\$", color="#377eb8" )
        ax3.plot( angle, M14/scale, label="\$M_{14}\$", color="#4daf4a" )
        ax3.plot( angle, M41/scale, label="\$M_{41}\$", color="#ff7f00" )
        ax3.legend( loc="best", frameon=False, labelspacing=0.05 )
        ax3.set_ylabel("\$\\times 10^{%d}\$"%(log10) )
        ax3.spines["right"].set_visible(False)
        ax3.spines["top"].set_visible(False)
        ax3.xaxis.set_ticks_position("bottom")
        ax3.yaxis.set_ticks_position("left")

        maxVals = np.max( [np.max(M23),np.max(M32),np.max(M24),np.max(M42)] )
        log10 = int( np.log10(maxVals) )
        if ( log10 < 0 ):
            log10 -= 1
        scale = 10**log10
        ax4.plot( angle, M23/scale, label="\$M_{23}\$", color="#e41a1c" )
        ax4.plot( angle, M32/scale, label="\$M_{32}\$", color="#377eb8" )
        ax4.plot( angle, M24/scale, label="\$M_{24}\$", color="#4daf4a" )
        ax4.plot( angle, M42/scale, label="\$M_{42}\$", color="#ff7f00" )
        ax4.set_ylabel("\$\\times 10^{%d}\$"%(log10) )
        ax4.legend( loc="best", frameon=False, labelspacing=0.05 )
        ax4.spines["right"].set_visible(False)
        ax4.spines["top"].set_visible(False)
        ax4.xaxis.set_ticks_position("bottom")
        ax4.yaxis.set_ticks_position("left")



def main( argv ):
    if ( len(argv) != 1 ):
        print ("Usage: python computeMueller.py hdf5file.h5")
        return 1

    mueller = MuellerMatrix()
    with h5.File(argv[0],'r') as hf:
        mueller.load(hf)

    mueller.plot()
    mueller.plotStokesVectorVSPhi( 3, -1 )
    mueller.degreeOfLinearPolarizationUnpolarized( [np.pi/4.0, 3.0*np.pi/4.0] )
    mueller.degreeOfLinPolPhiAveraged()
    mueller.muellerMatrixSingleElements( np.pi/4.0 )
    plt.show()

if __name__ == "__main__":
    main( sys.argv[1:] )
