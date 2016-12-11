import sys
#sys.path.append("../FresnelFDTD")
sys.path.append("../")
#import mplLaTeX as ml
import matplotlib as mpl
#mpl.rcParams.update(ml.params)
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 28
mpl.rcParams["axes.unicode_minus"]=False
import numpy as np
import json
import h5py as h5
from matplotlib import pyplot as plt
import colorScheme as cs
DELTA = 4.14E-5 # Salditt et al
BETA = 3.45E-6 # Salditt et al
import subprocess

def damping( u, v, k, width, eigenvalue, R ):
    q = np.sqrt( k**2 - eigenvalue )
    print q
    return np.exp(-BETA*v*k**2/q*(1.0+2.0*u/R))

def plot2DInRealCrd( u, transverseProfile, stat ):
    R = stat["outerRadius"]
    v = np.linspace(0.0, 5E5, 1001)
    U, V = np.meshgrid(u,v)
    A = damping(U, V, stat["wavenumber"], stat["width"], stat["solver"]["eigenvalue"], R)
    energy = np.zeros(U.shape)
    for i in range(0, U.shape[1]):
        if (( U[0,i] > 0.0 ) or ( U[0,i] < -stat["width"])):
            energy[:,i] = transverseProfile[i]*A[:,i]
        else:
            energy[:,i] = transverseProfile[i]

    plt.contourf(U,V/1000.0, energy, 100, cmap="gist_heat", norm=mpl.colors.LogNorm())
    plt.xlabel("$u$ (nm)")
    plt.ylabel("$v \; (\mathrm{\mu m})$")
    fname = "Figures/profile2D_uvplane.jpeg"
    plt.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))

    plt.clf()

    XYcompl = R*np.exp((U+1j*V)/R)
    transverse = XYcompl.real
    longitudinal = XYcompl.imag

    plt.contourf(longitudinal/1000.0,transverse,energy, 100, cmap="gist_heat", norm=mpl.colors.LogNorm())
    plt.gca().set_axis_bgcolor("#3E0000")
    plt.xlabel("$z \; (\mathrm{\mu m}$)")
    plt.ylabel("$x$ (nm)")
    fname = "Figures/profile2D_xyplane.jpeg"
    plt.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))
    return longitudinal, transverse, energy

def computeEffectiveFieldAbsorption( x, profile, wgstart, wgend ):
    totalpower = np.sum(profile**2)
    indxstart = np.argmin( np.abs(x-wgstart) )
    indxend = np.argmin( np.abs(x-wgend) )
    inside = np.sum(profile[indxstart:indxend]**2)
    outside = totalpower-inside
    return BETA*outside/totalpower

def addShadedBkg(ax, wgstart, wgend):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    width1 = wgstart-xmin
    height = ymax-ymin
    width2 = xmax-wgend
    color = "#d3d3d3"
    R1 = mpl.patches.Rectangle((xmin,ymin), width1, height, facecolor=color, edgecolor="none")
    R2 = mpl.patches.Rectangle((wgend,ymin), width2, height, facecolor=color, edgecolor="none")
    ax.add_patch(R1)
    ax.add_patch(R2)
    return ax

def main(argv):
    fname = ""
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        else:
            print ("Unknown argument %s"%(fname))
            return

    if ( fname == "" ):
        print ("No file specified")
        return

    infile = open( fname , "r" )
    stat = json.load(infile)
    infile.close()

    modeStart = 0
    modeEnd = 3
    nModes = 20
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    mxColors = len(cs.COLORS)
    absorption = np.zeros(nModes)
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3"]
    signs = "+---" # TODO: Make this more generic. Input from terminal?
    with h5.File(stat["solutionfile"], "r") as hf:
        for i in range(0, nModes):
            data = np.array( hf.get("mode%d"%(i)) )
            u = np.linspace( stat["solver"]["xmin"], stat["solver"]["xmax"], len(data))
            absorption[i] = computeEffectiveFieldAbsorption( u, data, -stat["width"], 0.0)
            print "Absorption mode %d: %.2E"%(i+1,absorption[i])
            if ( i >= modeStart ) and ( i < modeEnd ):
                if ( signs[i] == "+" ):
                    pre = 1.0
                else:
                    pre = -1.0
                ax.plot( u, pre*data, color=colors[i%4], label="%d"%(i+1))

    ymin, ymax = ax.get_ylim()
    ax = addShadedBkg( ax, -stat["width"], 0.0)
    ax.set_xlabel("$u$ (nm)")
    ax.set_ylabel("Intensity (a.u.)")
    ax.legend(loc="upper right", frameon=False, labelspacing=0.05)
    figname = "Figures/profile.svg"
    psname = "Figures/profile.ps"
    fig.savefig(figname, bbox_inches="tight")
    subprocess.call(["inkscape", "--export-ps=%s"%(psname), "--export-latex", figname])
    print ("Figure written to %s"%(figname))
    plt.show()

    # Plot potential
    with h5.File(stat["potentialFname"], "r") as hf:
        potential = np.array( hf.get(stat["potentialLabel"]))

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    xPot = np.linspace( stat["potentialXmin"], stat["potentialXmax"], len(potential))
    ax.plot(xPot, potential, 'k')
    ax.set_xlabel("$u$ (nm)")
    ax.set_ylabel("Potential (nm$^{-2}$)")
    fname = "Figures/potential.pdf"
    fig.savefig( fname, bbox_inches="tight")
    print ("Figure written to %s"%(fname))

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(absorption, ls="steps", color="black")
    ax.set_xlabel("Mode number")
    ax.set_ylabel("Effective absoroption coefficient")
    fname = "Figures/absorption.pdf"
    fig.savefig(fname, bbox_inches="tight")
    print ("Absorption written to %s"%(fname))

    # TODO: Adapt these functions to the new json files
    #longitudinal, transverse, energy = plot2DInRealCrd(u, data, stat )
if __name__ == "__main__":
    main(sys.argv[1:])
