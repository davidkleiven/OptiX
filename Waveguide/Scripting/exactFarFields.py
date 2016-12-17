import numpy as np
from scipy import optimize

class ExactFarFieldDefault:
    def __init__(self):
        self.name = "Default farfield"
        self.color = "#053061"

    def plot( self, ax, qmin, qmax ):
        return

    def fit( self, x, y ):
        return

    def normalize( self, data ):
        return

class TwoMirrorGaussians(ExactFarFieldDefault):
    def __init__(self):
        ExactFarFieldDefault.__init__(self)
        self.angle = 0.0
        self.wavenumber = 0.0
        self.wgwidth = 100.0
        self.name = "TwoMirrorGaussians"
        self.nmodes = 4
        self.amplitudes = np.array([1.0,1.0,1.0,1.0])

    def hermitteSeries( self, freq ):
        amp = self.amplitudes + 1j*np.zeros(len(self.amplitudes))
        for i in range(0,len(amp)):
            amp[i] *= (-1j)**i
        return np.polynomial.hermite.hermval(  freq*self.wgwidth, amp )

    def gaussHermitte( self, freq ):
        return np.exp(-0.5*self.wgwidth**2 *freq**2)*self.hermitteSeries( freq )

    def farField( self, freq ):
        kx = self.wavenumber*np.sin(self.angle)
        ff = np.zeros(len(freq)) + 1j*np.zeros( len(freq) )
        ff = ( self.gaussHermitte(freq+kx) + self.gaussHermitte(freq-kx) )
        return np.abs(ff)**2

    def fitFunc( self, freq, width, amp0, amp1, amp2, amp3 ):
        self.wgwidth = width
        self.amplitudes = np.array( [amp0, amp1, amp2, amp3] )
        return self.farField( freq )

    def fit( self, x, y ):
        p0 = np.zeros( len(self.amplitudes) +1 )
        p0[0] = self.wgwidth
        p0[1:] = self.amplitudes
        optimize.curve_fit( self.fitFunc, x, y, p0=p0 )
        print ("Fitted with: %.2E nm"%(self.wgwidth))

    def plot( self, ax, qmin, qmax ):
        # Initialize x-axis values
        xmin, xmax = ax.get_xlim()

        q = np.linspace( qmin, qmax, 1000 )
        x = np.linspace( xmin, xmax, len(q))
        ft = self.farField( q )
        ax.plot( x, ft, color=self.color, lw=3, ls="--", label="Gaussian" )
        return ax

class YoungSlit(ExactFarFieldDefault):
    def __init__(self):
        ExactFarFieldDefault.__init__(self)
        self.width = 100.0
        self.separation = 100.0
        self.dataMax = 1.0

    def farField( self, freq ):
        return 2.0*self.width*np.sinc(self.width*freq*0.5/np.pi)*np.cos( self.separation*freq*0.5 )

    def plot( self, ax, qmin, qmax ):
        xmin, xmax = ax.get_xlim()

        q = np.linspace( qmin, qmax, 1000 )
        x = np.linspace( xmin, xmax, len(q))
        ft = self.farField( q )**2
        ft *= self.dataMax/np.max(ft)

        # Normalize
        ax.plot( x, ft, color=self.color, lw=1, label="Slit" )
        return ax

    def normalize( self, data ):
        self.dataMax = np.max(data)
