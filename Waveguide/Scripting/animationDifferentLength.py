from __future__ import print_function
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import h5py as h5
import subprocess
import sys

try:
    import colormaps as cmaps
    cmap = cmaps.viridis
except:
    cmap = "viridis"

def main( argv ):
    fname = ""
    folder = ""
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        elif ( arg.find("--folder=") != -1 ):
            folder = arg.split("--folder=")[1]

    if ( fname == "" ):
        print ("No filename specified!")
        return
    elif ( folder == "" ):
        print ("No folder specified!")
        return

    with h5.File(fname, 'r') as hf:
        framenr = 0
        for i in range(0, len(hf.keys())):
            key = "intensity%d"%(i)
            dset = hf.get(key)
            length = float(dset.attrs["position"])
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            intensity = np.array(dset)
            thetaMin = float(dset.attrs["thetaMin"])
            thetaMax = float(dset.attrs["thetaMax"])
            phiMin = float(dset.attrs["phiMin"])
            phiMax = float(dset.attrs["phiMax"])
            extent = [phiMin, phiMax, thetaMin, thetaMax]
            im = ax.imshow(intensity, extent=extent, origin="lower", norm=mpl.colors.LogNorm(), cmap=cmap)
            ax.set_aspect( (extent[1]-extent[0])/(extent[3]-extent[2]))
            ax.set_xlabel("Exit angle (deg)")
            ax.set_ylabel("Incident angle (deg)")
            ax.set_title("Waveguide length: %.5f mm"%(length/1E6))
            figname = folder+"/frame%d.png"%(framenr+10)
            fig.savefig(figname)
            print("Figure written to %s"%(figname), end='\r')
            sys.stdout.flush()
            framenr += 1
            plt.close()
        print ("\n")

    #subprocess.call(["avconv", "-r", "1", "-i", folder+"/frame%%d.png", "-r","30", folder+"/anim.mp4"])

if __name__ == "__main__":
    main(sys.argv[1:])
