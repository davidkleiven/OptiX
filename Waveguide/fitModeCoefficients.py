import sys
sys.path.append("../FresnelFDTD")
sys.path.append("../")
import mplLaTeX as ml
import matplotlib as mpl
mpl.rcParams.update(ml.params)
import numpy as np
import json
import h5py as h5
from matplotlib import pyplot as plt
from scipy import stats
import colorScheme as cs

def main(argv):
    fname = ""
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        else:
            print ("Unknown argument %s"%(arg))

    with h5.File(fname, 'r') as hf:
        dset = hf.get(hf.keys()[0])
        coeff = np.array ( dset )
        zmin = dset.attrs.get("z0")
        zmax = dset.attrs.get("z1")
        fduid = int( dset.attrs.get("FDuid") )
        eiguid = int( dset.attrs.get("EIGuid") )

    ncoeff = coeff.shape[0]
    print ncoeff
    z = np.linspace( zmin, zmax, coeff.shape[1] )/1E6
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ncoeff = 6
    for i in range(0, ncoeff):
        currentCoeff = coeff[i,:]
        fitStart = int( len(currentCoeff)/2 )
        fitEnd = int( 3*len(currentCoeff)/4)
        slope, interscept, pvalue, rvalue, stderr = stats.linregress( z[fitStart:fitEnd], \
        np.log(np.abs(currentCoeff[fitStart:fitEnd])))
        fit = np.exp(interscept)*np.exp(slope*z)
        ax.plot( z, np.abs(currentCoeff), marker=",", color="black")
        ax.plot( z, fit, color=cs.COLORS[i], label="%d"%(i+1))
        print "Damping length: %.2E mm"%(-1.0/slope)

    ax.set_xlabel("$z$ (mm)")
    ax.set_ylabel("Mode coefficient")
    ax.set_yscale("log")
    ax.legend(loc="upper right", frameon=False)
    fname = "Figures/decayModeCoeff_fd%d_eig%d.jpeg"%(fduid, eiguid)
    fig.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))

if __name__ == "__main__":
    main( sys.argv[1:] )
