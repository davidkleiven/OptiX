import numpy as np
import circularWG as cwg

class CurvedWGCartesian:
    def __self__(self):
        cwg.CurvedWaveguide.__init__(self)
        self.x = []
        self.z = []
        self.R = 1.0
        self.width = 1.0

    def isInsideGuide( self, x, z ):
        isAbove = np.sqrt( (self.R+x)**2 + z**2 ) > self.R
        isBelow = (self.R+x)**2 + z**2 < (self.R + self.width)**2
        return isAbove and isBelow

    def normalVector( self, x, z):
        theta = np.arctan(z/x)
        nhat = np.zeros(2)
        nhat[0] = np.cos(theta)
        nhat[1] = np.sin(theta)
        return nhat

    def reflect( self, x, z, direction ):
        normalComponent = np.dot( direction, self.normalVector(x,z) )
        tangComponent = direction - normalComponent
        return tangComponent - normalComponent

    def solve( self, x0, z0, dir, xmin, xmax, zmin, zmax, N ):
        self.x.append(x0)
        self.z.append(z0)
        z = z0
        x = x0
        dx = (xmax-xmin)/N
        dz = (zmax-zmin)/N
        while ( z < zmax and x < xmax and x > xmin and z > zmin):
            x += direction[0]*dx
            z += direction[1]*dz
            if ( not isInsideGuide(x,z) ):
                self.x.append(x)
                self.z.append(z)
                direction = self.reflect( x, z, direction)

    def plot( self, ax ):
        ax.plot( self.x, self.z, color="white")
