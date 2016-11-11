import numpy as np
from matplotlib import pyplot as plt
import h5py as h5

class Mode:
    def __init__(self):
        self.profile = []
        self.eigenvalue = []
        self.xmin = 0.0
        self.xmax = 0.0

class Eigenmodes:
    def __init__( self ):
        self.modes = []
        self.beta = None
        self.delta = None
        self.width = None

    def read( self, h5file ):
        indx = 0
        while ( True ):
            dsetname = "mode%d"%(indx)
            if ( not dsetname in h5file.keys() ):
                print ("Read %d datasets from the h5 file"%(indx))
                if ( self.beta is None ):
                    raise( "Could not read beta from h5file!")
                elif ( self.delta is None ):
                    raise( "Could not read delta from h5file!")
                elif ( self.width is None ):
                    raise( "Could not read the waveguide width from the h5file!")
                return

            mode = Mode()
            dset = h5file.get( dsetname )
            mode.profile = np.array( dset )
            eigval = dset.attrs.get( "eigenvalue" )
            xmin = dset.attrs.get( "xmin" )
            xmax = dset.attrs.get( "xmax" )

            if ( self.beta is None ):
                self.beta = dset.attrs.get( "beta" )
            if ( self.delta is None ):
                self.delta = dset.attrs.get( "delta" )
            if ( self.width is None ):
                self.width = dset.attrs.get( "width" )

            if ( eigval is None ):
                raise ("No eigenvalue attached to the file")
            elif ( xmin is None ):
                raise ( "No attribute named xmin for dset %s!"%(dsetname))
            elif ( xmax is None ):
                raise ( "No attribute named xmax for dset %s!"%(dsetname))
            mode.eigenvalue = float( eigval )
            mode.xmin = float(  xmin )
            mode.xmax = float( xmax )
            self.modes.append( mode )
            indx += 1

    def integrateMode( self, modenumber, x0, xmax ):
        N = len( self.modes[modenumber].profile )
        dx = ( self.modes[modenumber].xmax - self.modes[modenumber].xmin )/N

        # Locate the start indices
        x = np.linspace( self.modes[modenumber].xmin, self.modes[modenumber].xmax, N)
        xstart = np.argmin( np.abs( x-x0) )
        xEnd = np.argmin( np.abs( x-xmax) )

        if ( xstart < 0 ):
            xstart = 0
        if ( xEnd > N ):
            xEnd = N

        return dx*np.sum( self.modes[modenumber].profile[xstart:xEnd]**2 )

    def effectiveAbsorption( self ):
        effAbs = np.zeros( len( self.modes) )
        print ("Note: Assuming that the waveguide starts at x=0 and ends at x=width!")
        for i in range( 0, len(effAbs) ):
            total = self.integrateMode( i, -2*self.width, 4.0*self.width )
            inside = self.integrateMode( i, 0.0, self.width )
            outside = total - inside
            effAbs[i] = self.beta*outside/total
        return effAbs

    def propagationConstants( self, k0 ):
        # Subtract off the "main" propagation constant k
        # Prop const[n] = beta[n]-k_0 = -0.5*E/k_0^2
        propConst = np.zeros( len( self.modes ) )
        for i in range(0, len(self.modes) ):
            propConst[i] = -0.5*self.modes[i].eigenvalue/k0**2
        return propConst

    def computeInitialCoefficient( self, amplitudeIn ):
        # Amplitude in is a funciton with one argument x
        coeff = np.zeros( len(self.modes) )
        for i in range(0, len(self.modes) ):
            mode = self.modes[i]
            x = np.linspace( mode.xmin, mode.xmax, len(mode.profile) )
            dx = x[1]-x[0]
            coeff[i] = np.sum( mode.profile*amplitudeIn(x) )*dx
        return coeff

    def contour( self, coeff, propConst, absorption, wglength ):
        Nx = len( self.modes[0].profile )
        Nz = 2000
        z = np.linspace( 0.0, wglength, Nz )
        cmap = "viridis"
        field = np.zeros( (Nx,Nz) )+1j*np.zeros( (Nx,Nz) )
        #for i in range( 0, len(self.modes) ):
        for i in range( 0, 13 ):
            decayPart =  coeff[i]*np.exp(1j*propConst[i]*z)*np.exp(-absorption[i]*z)
            print decayPart.shape
            field += np.outer( self.modes[i].profile, decayPart )
        field = np.abs( field )

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        im = ax.imshow( field )
        plt.show()
