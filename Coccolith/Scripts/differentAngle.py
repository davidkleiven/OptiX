import sys
sys.path.append("Scripts")
import numpy as np
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["axes.unicode_minus"] = False
mpl.rcParams["font.size"] = 18
import plotSpectrum as ps
from matplotlib import pyplot as plt

def main():
    fnames = ["data/CaCO3CrossSectionLong.h5", "data/CaCO3Coccolith5degY.h5", "data/CaCO3Coccolith10degY.h5", "data/CaCO3Coccolith15deg.h5"]
    angles = [0,5,10,15]
    specObjs = []
    maxWavelength = 800.0

    for fname in fnames:
        specObjs.append( ps.initSpectrum(fname) )

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3"]
    for i in range(0, len(fnames)):
        crossSection = specObjs[i].getScatteringCrossSecion()
        f = np.linspace( specObjs[i].freqmin, specObjs[i].freqmax(), len(crossSection) )
        wavelength = specObjs[i].voxelSize/f
        crossSection = crossSection[ wavelength < maxWavelength ]
        wavelength = wavelength[ wavelength < maxWavelength ]
        ax.plot( wavelength, crossSection, color=colors[i], label="\$\SI{%d}{\degree}\$"%(angles[i]))
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Scattering cross section (\$\SI{}{\micro\meter}\$)")

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.legend( loc="best", labelspacing=0.05, frameon=False )
    plt.show()
if __name__ == "__main__":
    main()
