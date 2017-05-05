import sys
import farFieldExactSphere as ffes
import numpy as np
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 28
mpl.rcParams["axes.unicode_minus"] = False
from matplotlib import pyplot as plt
import h5py as h5

def main( argv ):
    if ( len(argv) != 1 ):
        print ("Usage: python multiLayeredSphere.py hdffile.h5")
        return 1
    layers = ffes.LayeredSphere()
    layers.radii = np.array([1541.0,1597.0,1630.0,1691.0])
    layers.delta = np.array([5.46e-6,3.39e-5,6.27e-5,1.48e-5])

    # Read data
    with h5.File( argv[0], 'r' ) as hf:
        dset = hf.get("data/farField")
        qmin = dset.attrs["qmin"]
        qmax = dset.attrs["qmax"]
        ff = np.array( dset )

    ffLine = ff[int(ff.shape[0]/2),:]
    q = np.linspace( qmin, qmax, len(ffLine) )
    formfactor = layers.formFactor( q )**2

    # Normalize the data
    ffLine *= np.sum(formfactor)/np.sum(ffLine)

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot( q, ffLine, color="#e41a1c", label="Num")
    ax.plot( q, formfactor, color="#377eb8", label="\$|F(q)|^2\$" )
    ax.set_yscale("log")
    ax.set_xlabel("\$q_x (\SI{}{\per\\nano\meter})\$")
    ax.set_ylabel("Intensity (a.u.)")
    ax.legend(loc="upper left", frameon=False, labelspacing=0.05)
    plt.show()

if __name__ == "__main__":
    main( sys.argv[1:] )
