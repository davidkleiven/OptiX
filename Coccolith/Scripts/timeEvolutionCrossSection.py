import numpy as np
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["axes.unicode_minus"] = False
mpl.rcParams["font.size"] = 18
from matplotlib import pyplot as plt

def main():
    fnameRef = "data/referencePlaneTimeEvolution.csv"
    fnameBox = "data/fluxBoxTimeEvolution.csv"

    ref = np.loadtxt( fnameRef, delimiter="," )
    box = np.loadtxt( fnameBox, delimiter="," )

    fig = plt.figure()
    ax = fig.add_subplot(1,2,1)
    for i in range(0, ref.shape[1] ):
        ax.plot( np.abs(box[:,i]/ref[:,i]) )
    ax.legend()

    N = 20
    ax2 = fig.add_subplot(1,2,2)
    delta = int( box.shape[0]/N )
    time = np.arange(0,box.shape[1])*500
    for i in range(0,N):
        indx = i*delta
        data = np.abs( box[indx,:]/ref[indx,:] )
        ax2.plot( time, data )
        ax2.set_xlabel("Time")
    plt.show()

if __name__ == "__main__":
    main()
