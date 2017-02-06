import sys
import numpy as np
import matplotlib as mpl
import h5py as h5
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 36
mpl.rcParams["axes.unicode_minus"]=False
from matplotlib import pyplot as plt

# Transitions:
# 1) z=0.69 mm
# 2) z=2.09 mm

# Call: python3 connectionConsitency.py datafile.h5
def main( argv ):
    #fname = "data/multipleWG219886.h5"
    fname = "data/multipleWG323693.h5"
    fname = argv[0]
    with h5.File(fname,'r') as hf:
        group = hf.get("/data")
        xmin = group.attrs.get("xmin")
        xmax = group.attrs.get("xmax")
        zmin = group.attrs.get("zmin")/1E6
        zmax = group.attrs.get("zmax")/1E6
        amplitude = np.array( hf.get("/data/amplitude") )

    field = amplitude.T
    x = [-75.0,-50.0,-25.0 ]
    N = len( field[:,0] )
    indx = (x-xmin)*N/(xmax-xmin)
    z = np.linspace( zmin,zmax,len(field[0,:]) )

    colors = ["#66c2a5", "#fc8d62", "#8da0cb"]
    colors = ["#e41a1c", "#377eb8", "#4daf4a"]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    counter = 0
    for i in indx:
        ax.plot(z, np.real( field[i,:] ), color=colors[counter], label=x[counter] )
        counter += 1
    ax.spines["right"].set_visible( False )
    ax.spines["top"].set_visible( False )
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    #ax.legend(loc="upper right", labelspacing=0.01, frameon=False)
    plt.show()

if __name__ == "__main__":
    main( sys.argv[1:])
