import numpy as np
import skfmm
from matplotlib import pyplot as plt

class Lens:
    def __init__( self ):
        self.R1 = 100.0
        self.R2 = 100.0
        self.centerX = 0.0
        self.thickness = 40.0
        self.refractiveIndex = 1.5

    def isInside( self, x, y ):
        centerFirst = self.R1 - self.thickness/2.0
        centerSecond = -self.R2 + self.thickness/2.0
        r1 = np.sqrt( (x-centerFirst)**2 + y**2 )
        r2 = np.sqrt( (x-centerSecond)**2 + y**2 )
        return ( r1 < self.R1 ) and ( r2 < self.R2 )

    def plotRefractiveIndex( self, x, y ):
        refr = np.zeros((len(x),len(y)) )
        for i in range(0,len(x)):
            for j in range(0, len(y) ):
                refr[i,j] = self.isInside( x[j], y[i] )
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.imshow( refr )
        plt.show()

def main():
    n = 1.333
    N = 1000
    R = 1000
    xmin = -100.0
    xmax = 100.0
    x = 0.7*np.linspace(xmin,xmax,N)
    y = np.linspace(xmin,xmax,N)
    lens = Lens()
    lens.plotRefractiveIndex(x,y)
    speed = np.ones((N,N))
    for i in range(0,N):
        for j in range(0,N):
            if ( lens.isInside(x[i],y[j]) ):
                speed[i,j] = 1.0/lens.refractiveIndex

    # Define zero contour
    phi = np.ones(speed.shape)
    phi[1,1:phi.shape[1]-1] = -1
    travelTime = skfmm.travel_time(phi, speed)

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    X, Y = np.meshgrid(x,y)
    im = ax.contour(X,Y,travelTime, N=100)
    im = ax.imshow(travelTime)
    grad = np.gradient(travelTime)
    maxQuivers = 20
    X = X[::int(X.shape[0]/maxQuivers),::int(X.shape[1]/maxQuivers)]
    Y = Y[::int(Y.shape[0]/maxQuivers),::int(Y.shape[1]/maxQuivers)]
    #g1 = grad[1]
    #g1 = g1[::int(g1.shape[0]/maxQuivers),::int(g1.shape[1]/maxQuivers)]
    #g0 = grad[0]
    #g0 = g0[::int(g0.shape[0]/maxQuivers),::int(g0.shape[1]/maxQuivers)]

    #ax.quiver(X,Y,g1,g0)

    fig.colorbar(im)
    plt.show()

if __name__ == "__main__":
    main()
