import sys
sys.path.append("../FresnelFDTD")
sys.path.append("../")
import colorScheme as cs
import mplLaTeX as ml
import matplotlib as mpl
#mpl.rcParams.update(ml.params)
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["axes.linewidth"] = 0.1
mpl.rcParams["font.size"] = 28
import numpy as np
import h5py as h5
from matplotlib import pyplot as plt
import subprocess

FNAMES = ["data/singleCurvedWG133595_trans.h5", "data/singleCurvedWG558175_trans.h5", "data/transmissionExport.h5"]
LABELS = ["Cylindrical", "Cartesian", "Eigenmodes"]
DSET_NAMES = ["transmission", "transmission", "transmission"]
FIGNAME = "Figures/compareTransmission"
PSNAME = FIGNAME+".ps"
FIGNAME += ".svg"

def main():
    colors = ["#e41a1c", "#377eb8", "#4daf4a"]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(0, len(FNAMES)):
        print (FNAMES[i])
        with h5.File(FNAMES[i], 'r') as hf:
            dset = hf.get(DSET_NAMES[i])
            data = np.array( dset )
            zmin = dset.attrs.get("zmin")
            zmax = dset.attrs.get("zmax")

        z = np.linspace(zmin,zmax, len(data))
        ax.plot( z/1E6, np.log(data), color=colors[i], label=LABELS[i])
    ax.set_xlabel("\$z\$ (mm)")
    ax.set_ylabel("\$\ln T\$")
    ax.legend(loc="upper right", frameon=False)
    fig.savefig(FIGNAME, bbox_inches="tight")
    print ("Figure written to %s"%(FIGNAME))
    print ("Exporting to ps")
    subprocess.call(["inkscape", "--export-ps=%s"%(PSNAME), "--export-latex", FIGNAME])
    print ("PS-files written to %s"%(PSNAME))

if __name__ == "__main__":
    main()
