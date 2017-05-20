import sys
from scipy import io
import numpy as np
from matplotlib import pyplot as plt
import h5py as h5

def main( argv ):
    if ( len(argv) != 1 ):
        print ("Usage: python compareScatterWithExp.py hdf5file.h5" )

    matfile = io.loadmat("data/experimental_scattering.mat")
    expI = np.array( matfile["I1"][0] )
    expQ = np.array( matfile["Q"][0] )/1E9
    expI = expI[:len(expQ)]

    h5fname = argv[0]
    with h5.File(h5fname, 'r' ) as hf:
        ffdset = hf.get("data/farField")
        qmin = ffdset.attrs["qmin"]
        qmax = ffdset.attrs["qmax"]
        farfield = np.array(ffdset[int(ffdset.shape[0]/2),:])

    q = np.linspace(qmin,qmax, len(farfield))
    farfield = farfield[q>0.0]
    q = q[q>0.0]
    q = q[q<np.max(expQ)]
    farfield = farfield[:len(q)]

    # Normalize
    sumFF = np.sum(farfield[int(len(farfield)/5):])
    sumExpI = np.sum(  expI[int(len(expI)/5):])
    farfield *= np.abs(sumExpI)*len(farfield)/( np.abs(sumFF)*len(expI))

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(expQ, expI, label="Exp" )
    ax.plot(q, farfield, label="Num" )
    ax.set_yscale("log")
    plt.show()

if __name__ == "__main__":
    main( sys.argv[1:])
