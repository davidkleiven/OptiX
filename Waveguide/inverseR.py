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

FILES = ["singleCurvedWG265156_trans.h5","singleCurvedWG963983_trans.h5","singleCurvedWG681603_trans.h5",
"singleCurvedWG827383_trans.h5", "singleCurvedWG116355_trans.h5", "singleCurvedWG507568_trans.h5",
"singleCurvedWG680743_trans.h5", "singleCurvedWG326866_trans.h5"]
R = np.array( [10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 150.0] )
FOLDER = "data"
zmin = 0.0
zmax = 400.0
POS = np.array( [100, 200, 300, 400.0] )
#POS = np.array( [200.0] )
def main():
    assert( len(R) == len(FILES) )
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3"]
    assert ( len(POS) == len(colors) )
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(0, len(POS)):
        transmission = np.zeros(len(R))
        for j in range(0, len(FILES)):
            try:
                fname = FOLDER+"/"+FILES[j]
                with h5.File(fname, 'r') as hf:
                    t = np.array( hf.get("transmission"))
            except Exception as exc:
                print ("Error when reading HDF5 file %s"%(fname))
                print str(exc)
                return
            theta = POS[i]/(R[j]*1E3)
            theta_max = zmax/(R[j]*1E3)
            theta_min = zmin/(R[j]*1E3)
            indx = int( (theta-theta_min)*len(t)/(theta_max-theta_min) )
            if ( indx >= len(t) ):
                indx = -1
            transmission[j] = t[indx]
        invR = np.linspace(0.9*np.min(1/R), 1.05*np.max(1/R), 10)
        slope, interscept, pvalue, rvalue, stderr = stats.linregress(1.0/R, np.log(transmission))
        print ("Slope: %.2E mm"%(slope))
        print ("Interscept: %.2E"%(np.exp(interscept)))
        ax.plot( 1/R, np.log(transmission), marker="o", color="black", ms=7, ls="none", fillstyle="none")
        ax.plot( invR, interscept+slope*invR, color=colors[i], label="%d mm"%(POS[i]))
        ax.set_xlabel("\$R^{-1}\$ (mm\$^{-1}\$)")
        ax.set_ylabel("\$\ln T\$")
        ax.legend(bbox_to_anchor=(0.6,0.5), frameon=False, labelspacing=0.05)
        fname = "Figures/inverseRScaling.svg"
        fig.savefig(fname)
        psname = "Figures/inverseRScaling.ps"
        subprocess.call(["inkscape", "--export-ps=%s"%(psname), "--export-latex", fname])
        print ("Figure written to %s and %s"%(fname, psname))
    plt.show()

if __name__ == "__main__":
    main()
