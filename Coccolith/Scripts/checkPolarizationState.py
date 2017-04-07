import sys
import numpy as np
from matplotlib import pyplot as plt
import h5py as h5

def main( argv ):
    if ( len(argv) != 1 ):
        print ("Usage: python checkPolarizationState.py hdffile.h5")
        return 1

    with h5.File( argv[0], 'r' ) as hf:
        ey = np.array( hf.get("ey.r") )
        ez = np.array( hf.get("ez.r") )

    Nplots = 5
    Nx = 20
    Ny = 20
    for i in range(0,ey.shape[0], int(ey.shape[0]/Nplots) ):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        stepX = int(ey.shape[1]/Nx)
        stepY = int(ey.shape[2]/Ny)
        U = ey[i,::stepX,::stepY]
        V = ez[i,::stepX,::stepY]
        ax.imshow(U**2+V**2, cmap="inferno")
        ax.quiver(U,V,color="white")
        plt.show()

if __name__ == "__main__":
    main( sys.argv[1:] )
