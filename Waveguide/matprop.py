import numpy as np
import matplotlib as mpl
import subprocess
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 28
mpl.rcParams["axes.unicode_minus"]=False
from matplotlib import pyplot as plt
from scipy import stats

FNAME = "MatProp/indexRefrSiO2.txt"
FIGNAME = "Figures/refrIndexSiN3"
PSNAME = FIGNAME+".ps"
FIGNAME += ".svg"

def main():
    energy, delta, beta = np.loadtxt(FNAME, unpack=True, skiprows=2)
    energy /= 1000.0
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    every = 7
    fitStart = int( len(energy)/2 )
    slope, interscept, rvalue, pvalue, stderr = stats.linregress(np.log10(energy[fitStart:]), np.log10(delta[fitStart:]))
    fit = 10**(interscept)*energy**slope
    print ("Slope delta: %.1f"%(slope))

    ax.plot( energy[::every], delta[::every], color="black", marker='o', fillstyle="none", ls="none", label="\$\delta\$")
    ax.plot( energy, fit, color="black")
    ax.plot( energy[::every], beta[::every], color="black", marker="s", fillstyle="none", ls="none", lw=2, label="\$\\beta\$")
    slope, interscept, rvalue, pvalue, stderr = stats.linregress(np.log10(energy[fitStart:]), np.log10(beta[fitStart:]))
    fit = 10**(interscept)*energy**slope
    print ("Slope beta: %.1f"%(slope))
    ax.plot( energy, fit, color="black")
    ax.set_xlabel("Energy (\$\SI{}{\kilo\electronvolt}\$)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(right=10.0)
    ax.legend(loc="lower left", frameon=False, labelspacing=0.05)
    fig.savefig( FIGNAME )
    subprocess.call(["inkscape", "--export-ps=%s"%(PSNAME), "--export-latex", FIGNAME])
    print ("Figure written to %s"%(PSNAME))
    plt.show()

if __name__ == "__main__":
    main()
