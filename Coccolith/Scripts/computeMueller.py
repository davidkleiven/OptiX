import sys
sys.path.append("Scripts")
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import polarization as plz
import h5py as h5

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

    def load( self, hfile ):
        counter = 0
        key = "Run%d"%(counter)
        while (key in hfile.keys() ):
            pol = plz.Polarization()
            pol.readPhiAndTheta( hfile.get(key), 2 )
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
                im = axes[i*4+j].imshow( muellerAng[i][j].T[:,::-1], cmap=cmap, aspect="auto", vmin=minval, vmax=maxval)

        for i in range(0,16):
            axes[i].xaxis.set_ticks([])
            axes[i].yaxis.set_ticks([])

        cbar_ax = fig.add_axes([0.15, 0.95, 0.7, 0.05])
        fig.colorbar( im, orientation="horizontal", cax=cbar_ax )

    def plotStokesVectorVSPhi( self, runIndx, thetaIndx ):
        fig = plt.figure()
        ax1 = fig.add_subplot(1,2,1)
        ax1.plot( self.pol[runIndx].IphiTheta[thetaIndx,:], label="I")
        ax1.plot( self.pol[runIndx].QphiTheta[thetaIndx,:], label="Q")
        ax1.plot( self.pol[runIndx].UphiTheta[thetaIndx,:], label="U")
        ax1.plot( self.pol[runIndx].VphiTheta[thetaIndx,:], label="V")
        ax1.legend()

        ax2 = fig.add_subplot(1,2,2)
        ax2.plot( self.stokesInc[runIndx].Iinc, label="I")
        ax2.plot( self.stokesInc[runIndx].Qinc, label="Q")
        ax2.plot( self.stokesInc[runIndx].Uinc, label="U")
        ax2.plot( self.stokesInc[runIndx].Vinc, label="V")
        ax2.legend()
        plt.show()

def main( argv ):
    if ( len(argv) != 1 ):
        print ("Usage: python computeMueller.py hdf5file.h5")
        return 1

    mueller = MuellerMatrix()
    with h5.File(argv[0],'r') as hf:
        mueller.load(hf)

    mueller.plot()
    mueller.plotStokesVectorVSPhi( 3, -1 )
    print ( mueller.projectionMatrix(-1,-1) )
    plt.show()

if __name__ == "__main__":
    main( sys.argv[1:] )
