import numpy as np
import matplotlib as mpl
import subprocess
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 28
from matplotlib import pyplot as plt

FNAME_REAL = "MatPropData/refrIndexRealSiN3.csv"
FNAME_IMAG = "MatPropData/refrIndexImagSiN3.csv"
FIGNAME = "Figures/refrIndexSiN3"
PSNAME = FIGNAME+".ps"
FIGNAME += ".svg"

def main():
    energyRe, nre = np.loadtxt( FNAME_REAL, delimiter=",", unpack=True)
    energyIm, nim = np.loadtxt( FNAME_IMAG, delimiter=",", unpack=True)

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot( energyRe/1000.0, nre, color="black", label="\$\delta\$")
    ax.plot( energyIm/1000.0, nim, color="black", ls="--", lw=2, label="\$\\beta\$")
    ax.set_xlabel("Energy (\$\SI{}{\kilo\electronvolt}\$)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    #ax.legend(loc="upper right", frameon=False)
    ax.legend(bbox_to_anchor=(1.3,1), frameon=False)
    fig.savefig( FIGNAME )
    subprocess.call(["inkscape", "--export-ps=%s"%(PSNAME), "--export-latex", FIGNAME])
    print ("Figure written to %s"%(PSNAME))

if __name__ == "__main__":
    main()
