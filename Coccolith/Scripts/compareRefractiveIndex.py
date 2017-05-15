import numpy as np
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 18
mpl.rcParams["axes.unicode_minus"] = False
from matplotlib import pyplot as plt

def main():
    fnameWater = "Materials/water.csv"
    fnameCaCO3 = "Materials/CaCO3ne.csv"

    water = np.loadtxt( fnameWater, delimiter="," )
    caco3 = np.loadtxt( fnameCaCO3, delimiter="," )

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot( water[:,0]*1E3, water[:,1], label="\ce{H20}", color="#e41a1c" )
    ax.plot( caco3[:,0]*1E3, caco3[:,1], label="\ce{CaCO3}", color="#377eb8")
    ax.set_xlim(200,800)
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Refractive index")
    ax.legend(loc="best", frameon=False, labelspacing=True)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    plt.show()

if __name__ == "__main__":
    main()
