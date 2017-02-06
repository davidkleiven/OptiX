import sys
import numpy as np
from matplotlib import pyplot as plt
import h5py as h5

# Transitions:
# 1) z=0.69 mm
# 2) z=2.09 mm

def main():
    fname = "data/multipleWG219886.h5"
    with h5.File(fname,'r') as hf:
        group = hf.get("/data")
        xmin = group.attrs.get("xmin")
        xmax = group.attrs.get("xmax")
        zmin = group.attrs.get("zmin")/1E6
        zmax = group.attrs.get("zmax")/1E6
        amplitude = np.array( hf.get("/data/amplitude") )

    field = amplitude
    x = [-75.0,-50.0,-25.0 ]
    N = len( field[:,0] )
    indx = (x-xmin)*N/(xmax-xmin)
    z = np.linspace( zmin,zmax,len(field[0,:]) )

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in indx:
        ax.plot(z, np.real( field[i,:] ) )
    plt.show()

if __name__ == "__main__":
    main()
