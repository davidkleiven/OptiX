import sys
sys.path.append("../../FresnelFDTD")
import matplotlib as mpl
import mplLaTeX as ml
mpl.rcParams.update(ml.params)
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate

MSG = "Usage: python pattern.py --uid=<UID>\n"
FOLDER = "result" # UID will be appended
def main(argv):
    for arg in argv:
        if ( arg.find("--uid=") != -1 ):
            uid  = int(arg.split("--uid=")[1])
    fname = FOLDER+"%d"%(uid)+"/E_obs_scatt.txt"
    x, y, z, Exr, Exim, Eyr, Eyim, Ezr, Ezim = np.loadtxt(fname, unpack="True", skiprows=1)
    
    xInterp = np.linspace(np.min(x), np.max(x), 101)
    yInterp = np.linspace(np.min(x), np.max(x), 101)
    X, Y = np.meshgrid(xInterp, yInterp)
    
    intensity = Exr**2 + Exim**2 + Eyr**2 + Eyim**2 + Ezr**2 + Ezim**2
    intensity = interpolate.griddata((x,y), intensity,  (X,Y))

    plt.contourf( X,Y, intensity, 200,cmap="gist_heat")
    fname = "Fig/pattern2d%d.png"%(uid)
    plt.savefig(fname, bbox_inches="tight")
    print ("Figure written to %s"%(fname))
    
if __name__ == "__main__":
    main(sys.argv[1:])
