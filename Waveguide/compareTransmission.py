import sys
#sys.path.append("../FresnelFDTD")
sys.path.append("../")
import colorScheme as cs
#import mplLaTeX as ml
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 28
import numpy as np
import h5py as h5
import json
from matplotlib import pyplot as plt
from scipy import stats
import subprocess

def main( argv ):
    MSG = "Usage: python compareTransmission.py --file=<jsonfile> [--help]\n"
    MSG += "help: Print this message\n"
    MSG += "file: Json files describing the UIDs to plot"
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        elif ( arg.find("--help") != -1 ):
            print MSG
            return 0
        else:
            print ("Unknown argument %s"%(arg))
            return 0

    try:
        infile = open(fname,'r')
        param = json.load(infile)
        infile.close()
    except Exception as exc:
        print str(exc)
        print ("Error when opening/parsing file %s"%(fname))
        return 0

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    indx = 0
    markers = ["o", "o", "o"]
    fs = ["full", "full", "full"]
    color = ["#e41a1c", "#377eb8", "#4daf4a"]
    minOfMaxZ = np.inf
    for entry in param["entries"]:
        ctlfile = param["basename"]+"%d.json"%(entry["uid"])
        try:
            infile = open(ctlfile, 'r')
            stat = json.load(infile)
            infile.close()
        except Exception as exc:
            print str(exc)
            print ("Error when opening/parsing file %s"%(ctlfile))
            return 0

        with h5.File(stat["Transmission"]["file"], 'r') as hf:
            data = np.array( hf.get( "transmission" ) )

        try:
            crd = stat["waveguide"]["crd"]
        except:
            crd = "cartesian"

        # Reduce number of points in plot
        data = data[::param["numberOfPoints"]]
        z = np.linspace(stat["Transmission"]["zStart"], stat["Transmission"]["zEnd"], len(data))
        if ( crd == "cylindrical" ):
            z*= stat["waveguide"]["RadiusOfCurvature"]
        fitStart = np.argmin( np.abs(z/1E6 - 0.3) )
        zFit = z[fitStart:]
        dataFit = data[fitStart:]
        slope, interscept, rvalue, pvalue, stderr = stats.linregress(zFit,np.log(dataFit))
	zFit = np.linspace(0.4*np.max(z), 1.05*np.max(z), 11)

        print ("Damping length %s mm: %.2E mm"%(entry["label"],-1.0/(slope*1E6)))
        if ( stat["Transmission"]["zEnd"] < minOfMaxZ ):
            minOfMaxZ = stat["Transmission"]["zEnd"]
            ymin = np.min(np.log(data))
        ax.plot( z/1E6, np.log(data), marker=markers[indx], ms=7, color=color[indx], fillstyle=fs[indx], linestyle="None", label=entry["label"])
        ax.plot( zFit/1E6, interscept+slope*zFit, lw=0.5, color="black" )
        indx += 1

    ax.set_xlabel("\$z\$ (mm)")
    #ax.set_xlim(right=minOfMaxZ/1E6)
    ax.set_ylim(bottom=ymin-0.05)
    ax.set_ylabel("\$\ln Transmission \$")
    ax.legend(loc="upper right", frameon=False)
    fname = "Figures/"+param["figurename"]
    fig.savefig(fname, bbox_inches="tight")
    print ("Figure written to %s"%(fname))

    if ( fname[-3:] == "svg"):
        psname = fname[:-3]+"ps"
        subprocess.call(["inkscape", "--export-ps=%s"%(psname), "--export-latex", fname])
	print ("PS version written to %s"%(psname))

if __name__ == "__main__":
    main(sys.argv[1:])
