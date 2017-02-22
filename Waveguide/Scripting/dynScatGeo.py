import numpy as np

class Geometry:
    def __init__( self ):
        self.k = 40.0
        self.xmin = -1.0
        self.xmax = 1.0
        self.ymin = -1.0
        self.ymax = 1.0
        self.Nx = 128
        self.Ny = 128

    def getProjection( self, X, Y ):
        return 0.0

    def qlim( self ):
        qmax = 5.0*np.pi/(self.xmax-self.xmin)
        qmin = -qmax
        return qmin, qmax

    def farField( self ):
        x = np.linspace( self.xmin, self.xmax, self.Nx )
        y = np.linspace( self.ymin, self.ymax, self.Ny )
        X,Y = np.meshgrid(x,y)
        proj = self.getPhaseAccumulation( X, Y )
        del X, Y

        proj = np.fft.fft2( proj, s=(2048,2048) )
        proj = np.fft.fftshift( proj )
        return np.abs( proj )**2

class Sphere(Geometry):
    def __init__(self):
        Geometry.__init__(self)
        self.R = 500.0
        self.delta = 4.9E-5

    def getProjection( self, X, Y ):
        Rsq = X**2 + Y**2
        res = np.zeros( (Rsq.shape[0], Rsq.shape[1]) )
        res[Rsq<self.R*self.R] = self.delta*np.sqrt(self.R*self.R-Rsq[Rsq<self.R*self.R])
        return res

    def getPhaseAccumulation( self, X, Y ):
        Rsq = X**2 + Y**2
        res = np.zeros( (Rsq.shape[0], Rsq.shape[1]) )
        res[Rsq<self.R*self.R] = np.sqrt(self.R*self.R-Rsq[Rsq<self.R*self.R])
        return np.exp(-1j*self.k*self.delta*res ) - 1.0
