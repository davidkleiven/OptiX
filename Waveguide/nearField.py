import sys
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 28
import numpy as np
import h5py as h5
from matplotlib import pyplot as plt
import subprocess
try:
    import colormaps as cmaps
    cmap = cmaps.viridis
except:
    cmap = "viridis"

def main(argv):
    fname = ""
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        elif ( arg.find("--help") != -1 ):
            print ("Usage: python nearField.py --file=<hdf5data.h5> --help")
            print ("file: Filename")
            print ("help: Pinr this message")
        else:
            print ("Unknown argument %s"%(fname))
            return

    if ( fname == "" ):
        print ("No filename specified!")
        return

    with h5.File(fname, 'r') as hf:
        dset = hf.get("intensity")
        xmin = float( dset.attrs.get("xmin") )
        xmax = float( dset.attrs.get("xmax") )
        zmin = float( dset.attrs.get("zmin") )
        zmax = float( dset.attrs.get("zmax") )
        uid = int( dset.attrs.get("uid") )
        intensity = np.array(dset)

    fig = plt.figure()
    intensity /= np.max(intensity)
    print intensity.shape
    ax = fig.add_subplot(1,1,1)
    extent = [zmin/1000.0, zmax/1000.0, xmin, xmax]
    intensity = intensity.T
    im = ax.imshow( intensity, extent=extent, origin="lower", cmap=cmap )
    fig.colorbar(im)
    ax.set_aspect( np.abs( (extent[1]-extent[0])/(extent[3]-extent[2]) ))
    ax.set_xlabel("\$z (\SI{{\micro\meter}})\$")
    ax.set_ylabel("\$x\$ (nm)")
    ax.autoscale(False)
    fname = "Figures/nearField%d.svg"%(uid)
    psname = "Figures/nearField%d.ps"%(uid)
    jpname = "Figures/nearField%d.jpeg"%(uid)
    fig.savefig(fname)
    fig.savefig(jpname)
    subprocess.call( ["inkscape", "--export-ps=%s"%(psname), "--export-latex", fname] )
    print ("Figure written to %s, %s, %s"%(fname, psname, jpname))
    plt.show()

if __name__ == "__main__":
    main(sys.argv[1:])
