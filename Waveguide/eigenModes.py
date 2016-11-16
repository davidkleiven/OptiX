import sys
sys.path.append("../FresnelFDTD")
import mplLaTeX as mp
import matplotlib as mpl
#mpl.rcParams.update( mp.params )
mpl.rcParams['svg.fonttype'] = 'none'
import numpy as np
from matplotlib import pyplot as plt
import h5py as h5
try:
    import colormaps as cmaps
    plt.register_cmap(name="viridis", cmap=cmaps.viridis)
except:
    print ("Module colormaps not found")

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
        self.transmission = None
        self.nPropagatingModes = 0
        self.exportTransFname = ""
        self.radius = None
        self.normalizeContourPlots = True

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
                elif ( self.radius is None ):
                    raise( "Could not read the radius of curvature from the h5file!" )
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
            if ( self.radius is None ):
                self.radius = dset.attrs.get( "Rcurv" )

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

    def imaginaryPerturbation( self, modenumber ):
        N = len( self.modes[modenumber].profile )
        x = np.linspace( self.modes[modenumber].xmin, self.modes[modenumber].xmax, N)
        xstart = np.argmin( np.abs( x+self.width) )
        xEnd = np.argmin( np.abs(x) )
        total =  np.trapz( self.modes[modenumber].profile**2 )
        inside = np.trapz( self.modes[modenumber].profile[xstart:xEnd]**2 )
        outside = np.trapz( self.modes[modenumber].profile[:xstart]**2 ) + np.trapz( self.modes[modenumber].profile[xEnd:]**2)
        I1 = -2.0 *self.beta*(total-inside)/total
        I1 = -2.0 *self.beta*outside/total

        # Second term
        total2 = np.trapz( x*self.modes[modenumber].profile**2 )
        inside2 = np.trapz( x[xstart:xEnd]*self.modes[modenumber].profile[xstart:xEnd]**2 )
        I2 = -(4.0*self.beta /self.radius )*(total2-inside2)/total
        return I1+I2

    def integrateMode( self, modenumber, x0, xmax ):
        N = len( self.modes[modenumber].profile )
        dx = ( self.modes[modenumber].xmax - self.modes[modenumber].xmin )/N

        # Locate the start indices
        x = np.linspace( self.modes[modenumber].xmin, self.modes[modenumber].xmax, N)
        xstart = np.argmin( np.abs( x-x0) )
        xEnd = np.argmin( np.abs(x-xmax) )

        if ( xstart < 0 ):
            xstart = 0
        if ( xEnd > N ):
            xEnd = N
        return dx*np.trapz( self.modes[modenumber].profile[xstart:xEnd]**2 )

    def effectiveAbsorption( self ):
        k0 = 2.0*np.pi/0.1569 # Just for debuggin
        effAbs = np.zeros( len( self.modes) )
        print ("Note: Assuming that the waveguide starts at x=-width and ends at x=0.0!")
        for i in range( 0, len(effAbs) ):
            #total = self.integrateMode( i, -2*self.width, 4.0*self.width )
            #inside = self.integrateMode( i, -self.width, 0.0 )
            #outside = total - inside
            #effAbs[i] = self.beta*outside/total
            effAbs[i] =  -self.imaginaryPerturbation( i )/2.0 # Why 4.0?
        return effAbs

    def propagationConstants( self, k0 ):
        # Subtract off the "main" propagation constant k
        # Prop const[n] = beta[n]-k_0 = -0.5*E/k_0
        propConst = np.zeros( len( self.modes ) )
        for i in range(0, len(self.modes) ):
            propConst[i] = -0.5*self.modes[i].eigenvalue/k0
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

    def transmissionByIntegrateOverWG( self, coeff, propConst, absorption, k0, zmax ):
        z = np.linspace(0.0, zmax, 10001)
        T = np.zeros(len(z))
        Tfull = np.zeros(len(z))
        xmin = self.modes[0].xmin
        xmax = self.modes[0].xmax
        x = np.linspace(xmin, xmax, len(self.modes[0].profile))
        xstart = np.argmin( np.abs( x+self.width))
        xend = np.argmin( np.abs(x) )
        intensity = np.zeros((len(x),len(z))) + 1j*np.zeros((len(x), len(z)))
        for n in range(0, self.nPropagatingModes):
            intensity += self.fieldFromMode(n, coeff[n], propConst[n], absorption[n], k0, z)

        intensity = np.abs( intensity )**2
        for i in range(0, len(T)):
            T[i] = np.trapz( intensity[xstart:xend,i])
            Tfull[i] = np.trapz( intensity[:,i])

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        T /= T[0]
        Tfull /= Tfull[0]
        ax.plot( z/1E3, np.log(T), color="black" )
        ax.set_xlabel( "\$z (\mathrm{\mu m})\$" )
        ax.set_ylabel( "\$\ln T\$" )
        #fname = "Figures/transmissionModesIntegration.pdf"
        fname = "Figures/transmissionModesIntegration.svg"
        fig.savefig( fname, bbox_inches="tight")
        print ("Figure written to %s"%(fname))

        if ( self.exportTransFname != "" ):
            with h5.File( self.exportTransFname, 'w' ) as hf:
                dset1 = hf.create_dataset( "transmission", data=T)
                dset2 = hf.create_dataset( "transmissionFull", data=Tfull)
                dset1.attrs.create("zmin", 0.0)
                dset1.attrs.create("zmax", zmax)
                dset2.attrs.create("zmin", 0.0)
                dset2.attrs.create("zmax", zmax)
            print ("Transmission data written to %s"%(self.exportTransFname))

    def plotAbsorption( self, coeff, absCoeff, k0, zmax ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        z = np.linspace(0.0, zmax, 101 )
        T = np.zeros(len(z))
        #for i in range(0, len(coeff)):
        for i in range(0, self.nPropagatingModes):
            T += coeff[i]**2 *np.exp( -k0*z*absCoeff[i])

        T /= T[0]
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot( z/1E3, np.log(T), color="black" )
        ax.set_ylabel("$\ln T$")
        ax.set_xlabel("$z$ ($\mathrm{\mu m}$)")
        fname = "Figures/transmissionModes.pdf"
        fig.savefig(fname, bbox_inches="tight")
        print ("Figure written to %s"%(fname))

    def fieldFromMode( self, modenumber, coeff, prop, decay, k0, z ):
        beta = k0 + prop
        #decayPart =  coeff*np.exp(1j*prop*z)*np.exp(-0.5*decay*k0*z)
        decayPart =  coeff*np.exp(1j*prop*z)*np.exp(-decay*k0*z)
        return np.outer( self.modes[modenumber].profile, decayPart )

    def contour( self, coeff, propConst, absorption, k0, wglength ):
        cmap = "viridis"
        Nx = len( self.modes[0].profile )
        Nz = 2000
        z = np.linspace( 0.0, wglength, Nz )
        cmap = "viridis"
        field = np.zeros( (Nx,Nz) )+1j*np.zeros( (Nx,Nz) )
        for i in range( 0, self.nPropagatingModes ):
            field += self.fieldFromMode( i, coeff[i], propConst[i], absorption[i], k0, z)
        field = np.abs( field )**2

        appendName = ""
        if ( self.normalizeContourPlots ):
            field /= np.max(field)
            appendName = "normalized"

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        xmin = self.modes[0].xmin
        xmax = self.modes[0].xmax
        extent = [0.0, wglength/1E3, xmin, xmax]
        im = ax.imshow( field, extent=extent, aspect="equal", origin="lower", cmap=cmap)#, norm=mpl.colors.LogNorm() )
        fig.colorbar(im)
        ax.set_aspect( (extent[1]-extent[0])/(extent[3]-extent[2]) )
        ax.set_xlabel("$z (\mathrm{\mu m})$")
        ax.set_ylabel("$x$ (nm)")
        fname = "Figures/contour%s.jpeg"%(appendName)
        fig.savefig( fname, bbox_inches="tight", dpi=800)
        print ("Figure written to %s"%(fname))

        # Same on log scale
        fig = plt.figure()
        maxval = 0.1
        minval = 1E-6
        ax = fig.add_subplot(1,1,1)
        xmin = self.modes[0].xmin
        xmax = self.modes[0].xmax
        extent = [0.0, wglength/1E3, xmin, xmax]
        im = ax.imshow( field, extent=extent, aspect="equal", origin="lower", cmap=cmap, norm=mpl.colors.LogNorm(minval,maxval) )
        fig.colorbar(im)
        ax.set_aspect( (extent[1]-extent[0])/(extent[3]-extent[2]) )
        ax.set_xlabel("$z (\mathrm{\mu m})$")
        ax.set_ylabel("$x$ (nm)")
        fname = "Figures/contourLog%s.jpeg"%(appendName)
        fig.savefig( fname, bbox_inches="tight", dpi=800)
        print ("Figure written to %s"%(fname))
