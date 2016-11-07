import numpy as np
import circularWG as cwg
from matplotlib import pyplot as plt

class CurvedWGCartesian:
    def __init__(self):
        self.x = []
        self.z = []
        self.R = 1.0
        self.width = 1.0

    def isInsideGuide( self, x, z ):
        isAbove = np.sqrt( (self.R+x)**2 + z**2 ) > self.R
        isBelow = (self.R+x)**2 + z**2 < (self.R + self.width)**2
        return isAbove and isBelow

    def normalVector( self, x, z):
        theta = np.arctan(z/(self.R+x))
        nhat = np.zeros(2)
        nhat[0] = np.cos(theta)
        nhat[1] = np.sin(theta)
        return nhat

    def reflect( self, x, z, direction ):
        nhat = self.normalVector(x,z)
        normalComponent = np.dot( direction, nhat )
        tangComponent = direction - normalComponent*nhat
        direction = tangComponent - normalComponent*nhat
        return direction

    def solve( self, x0, z0, direction, xmin, xmax, zmin, zmax, N ):
        self.x.append(x0)
        self.z.append(z0)
        z = z0
        x = x0
        dt = (xmax-xmin)/N
        dz = (zmax-zmin)/N
        c = 1.0
        hasBeenInsideAfterReflcetion = True
        while ( z < zmax and x < xmax and x >= xmin and z >= zmin):
            xTrial = x + direction[0]*c*dt
            zTrial = z + direction[1]*c*dt
            if ( not self.isInsideGuide(xTrial,zTrial) ):
                self.x.append(x)
                self.z.append(z)
                direction = self.reflect( x, z, direction)

            x += direction[0]*c*dt
            z += direction[1]*c*dt

        self.x.append(x)
        self.z.append(z)

    def plot( self, ax ):
        self.z = np.array( self.z )
        self.x = np.array( self.x )
        ax.plot( self.z/1E3, self.x, color="white", lw=0.25)
