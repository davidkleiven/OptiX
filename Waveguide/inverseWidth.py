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
# "singleCurvedWG114123_trans.h5", w=10 nm, no modes
#FILES = ["singleCurvedWG301088_trans.h5","singleCurvedWG433090_trans.h5",
#"singleCurvedWG857974_trans.h5", "singleCurvedWG348856_trans.h5", "singleCurvedWG539750_trans.h5",
#"singleCurvedWG763629_trans.h5", "singleCurvedWG397296_trans.h5"]
#W = np.array( [20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 150.0] ) # in nm

FILES = ["singleCurvedWG326678_trans.h5","singleCurvedWG479615_trans.h5",
"singleCurvedWG943087_trans.h5", "singleCurvedWG963353_trans.h5", "singleCurvedWG524050_trans.h5",
"singleCurvedWG473836_trans.h5", "singleCurvedWG903590_trans.h5", "singleCurvedWG272447_trans.h5"]
W = np.array( [40.0, 60.0, 80.0, 100.0, 200.0, 400.0, 800.0, 1200.0] )
FOLDER = "data"
zmin = 0.0
zmax = 400.0
POS = np.array( [100, 200, 300, 400.0] )
R = 40.0
#POS = np.array( [200.0] )
def main():
    assert( len(W) == len(FILES) )

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(0, len(POS)):
        transmission = np.zeros(len(W))
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
        invWSq = np.linspace(0.9*np.min(1/W**2), 1.05*np.max(1/W**2), 10)

        slope, interscept, pvalue, rvalue, stderr = stats.linregress(1/W**2, np.log(transmission) )
        ax.plot( 1/W**2, np.log(transmission), marker="o", color="black", ms=7, ls="none", fillstyle="none")
        ax.plot( invWSq, interscept+slope*invWSq, color=cs.COLORS[i], label="%.1f mm"%(POS[i]/1E3))
        ax.set_xlabel("\$w^{-2}\$ (nm\$^{-2}\$)")
        ax.set_ylabel("Effective attenuation coefficient")
        ax.legend(bbox_to_anchor=(0.6,0.5), frameon=False, labelspacing=0.05)
        fname = "Figures/inverseWScaling.svg"
        fig.savefig(fname)
        psname = "Figures/invserseWScaling.ps"
        subprocess.call(["inkscape", "--export-ps=%s"%(psname), "--export-latex", fname])
        print ("Figure written to %s and %s"%(fname, psname))
    plt.show()

if __name__ == "__main__":
    main()
