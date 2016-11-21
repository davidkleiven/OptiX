import sys
sys.path.append("../FresnelFDTD")
sys.path.append("../")
import colorScheme as cs
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["axes.linewidth"] = 0.1
mpl.rcParams["font.size"] = 28
import numpy as np
import h5py as h5
from matplotlib import pyplot as plt
from scipy import stats
import subprocess

FILES = ["singleCurvedWG560421_trans.h5","singleCurvedWG964029_trans.h5",
"singleCurvedWG101998_trans.h5", "singleCurvedWG489293_trans.h5", "singleCurvedWG929243_trans.h5",
"singleCurvedWG881666_trans.h5", "singleCurvedWG139399_trans.h5","singleCurvedWG230470_trans.h5"]
DELTA = np.array( [1E-6, 3.162E-6, 5.623E-6, 1E-5, 1.778E-5, 3.162E-5, 5.6234E-5, 1E-4] )
FOLDER = "data"
zmin = 0.0
zmax = 400.0
POS = np.array( [100, 200, 300, 400.0] )
R = 40.0
#POS = np.array( [200.0] )
def main():
    assert( len(DELTA) == len(FILES) )

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(0, len(POS)):
        transmission = np.zeros(len(DELTA))
        for j in range(0, len(FILES)):
            try:
                fname = FOLDER+"/"+FILES[j]
                with h5.File(fname, 'r') as hf:
                    t = np.array( hf.get("transmission"))
            except Exception as exc:
                print ("Error when reading HDF5 file %s"%(fname))
                print str(exc)
                return
            theta = POS[i]/(R*1E3)
            theta_max = zmax/(R*1E3)
            theta_min = zmin/(R*1E3)
            indx = int( (theta-theta_min)*len(t)/(theta_max-theta_min) )
            if ( indx >= len(t) ):
                indx = -1
            transmission[j] = t[indx]
        invDelta = np.linspace(0.9*np.min(1/DELTA), 1.05*np.max(1/DELTA), 10)

        slope, interscept, pvalue, rvalue, stderr = stats.linregress(1/DELTA, np.log(transmission) )
        ax.plot( 1/DELTA, np.log(transmission), marker="o", color="black", ms=7, ls="none", fillstyle="none")
        ax.plot( invDelta, interscept+slope*invDelta, color=cs.COLORS[i], label="%.1f mm"%(POS[i]/1E3))
        ax.set_xlabel("\$\delta^{-1}\$")
        ax.set_ylabel("Effective attenuation coefficient")
        ax.legend(bbox_to_anchor=(0.6,0.5), frameon=False, labelspacing=0.05)
        fname = "Figures/inverseDeltaScaling.svg"
        fig.savefig(fname)
        psname = "Figures/invserseDeltaScaling.ps"
        subprocess.call(["inkscape", "--export-ps=%s"%(psname), "--export-latex", fname])
        print ("Figure written to %s and %s"%(fname, psname))
    plt.show()

if __name__ == "__main__":
    main()
