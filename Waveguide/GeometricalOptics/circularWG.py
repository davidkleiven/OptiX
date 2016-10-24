import sys
sys.path.append("../../FresnelFDTD")
import mplLaTeX as ml
import matplotlib as mpl
mpl.rcParams.update(ml.params)
import numpy as np
import matplotlib as mpl
from scipy import integrate
from matplotlib import pyplot as plt

class Ray:
    def __init__( self ):
        self.theta= None
        self.r = None

class CurvedWaveguide:
    def __init__(self):
        self.width = 0.0 # in nm
        self.R = 1.0 # in nm
        self.kappa = 1.0 # Assume that the effect of curvature is to change the refractive index to
        self.solution = []
        self.hist = None

    def rhs( self, rayvector, thetaDeg ):
        derivative = np.zeros(2)
        #print rayvector

        #if ( rayvector[0] > self.width ):
        #    rayvector[1] = -np.abs(rayvector[1])
        #elif ( rayvector[0] < 0 ):
        #    rayvector[1] = np.abs(rayvector[1])

        derivative[0] = rayvector[1]
        derivative[1] = self.kappa*self.R*(1.0 + self.force(rayvector[0]))
        return derivative

    def force( self, r ):
        factor = 100.0
        fwidth = 0.01
        return -factor/(1.0 + np.exp((self.width-r)/fwidth) )

    def solve( self, nRays, thetaMax, nTheta=1000 ):
        r0 = np.linspace(0.1, self.width-0.1, nRays )
        theta = np.linspace( 0.0, thetaMax, nTheta )
        for i in range(0, nRays):
            y = integrate.odeint( self.rhs, [r0[i],0.0], theta )
            solution = Ray()
            solution.theta = theta
            solution.r = y[:,0]
            self.solution.append( solution )

    def rayDensity( self ):
        print ("Computing ray density")
        for i in range(0, len(self.solution)):
            if ( i== 0):
                self.hist = np.histogram2d(self.solution[i].r, self.solution[i].theta, bins=100)[0]
            else:
                self.hist += np.histogram2d(self.solution[i].r, self.solution[i].theta, bins=100)[0]

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        minval = 0.0
        maxval = np.max(self.solution[0].theta)*180.0/np.pi
        extent = [0.0, maxval, 0.0,self.width]
        extend = [np.min(self.solution[0])]
        im = ax.imshow(self.hist, cmap="coolwarm", extent=extent, aspect=1, norm=mpl.colors.LogNorm())
        ax.set_xlabel("$\\theta$ (deg)")
        ax.set_ylabel("$r$ (nm)")
        ax.set_aspect( np.abs( (extent[1]-extent[0])/(extent[3]-extent[2]) ))
        fname = "Figures/rayDensity.jpeg"
        fig.savefig(fname, bbox_inches="tight", dpi=800)
        print ("Figure written to %s"%(fname))


    def plot( self ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        for i in range(0, len(self.solution) ):
            ax.plot(self.solution[i].theta*180.0/np.pi, self.solution[i].r, color="black", lw=0.2)

        ax.axhline(0.0, color="grey", lw=3)
        ax.axhline(self.width, color="grey", lw=3)
        ax.set_xlabel("$\\theta$ (deg)")
        ax.set_ylabel("$r$ (nm)")
        fname = "Figures/ray.pdf"
        fig.savefig(fname, bbox_inches="tight")
        print ("Figure written to %s"%(fname))
