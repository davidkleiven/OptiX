import plotSpectrum as psp
import numpy as np
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["axes.unicode_minus"] = False
mpl.rcParams["font.size"] = 18
from matplotlib import pyplot as plt

def main():
    #filledSideFile = "data/CaCO3Coccolith_20170505_1928.h5"
    filledSideFile = "data/CaCO3CoccolithEdgeFilled.h5"
    #filledCenterFile = "data/CaCO3Coccolith_20170505_1258.h5"
    filledCenterFile = "data/CaCO3Filled.h5"
    origFname = "data/CaCO3Coccolith_20170328_0804.h5"
    #origFname = "data/CaCO3CrossSectionLong.h5"

    spOrg = psp.initSpectrum( origFname )
    spSide = psp.initSpectrum( filledSideFile )
    spCenter = psp.initSpectrum( filledCenterFile )

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    scatCross = spOrg.getScatteringCrossSecion()
    f = np.linspace( spOrg.freqmin, spOrg.freqmax(), len( scatCross) )
    wavelength = spOrg.voxelSize/f
    ax.plot( wavelength, scatCross, label="Orig", color="#e41a1c")

    scatCross = spSide.getScatteringCrossSecion()
    f = np.linspace( spSide.freqmin, spSide.freqmax(), len( scatCross) )
    wavelength = spSide.voxelSize/f
    ax.plot( wavelength, scatCross, label="Side", color="#377eb8")


    scatCross = spCenter.getScatteringCrossSecion()
    f = np.linspace( spCenter.freqmin, spCenter.freqmax(), len( scatCross) )
    wavelength = spCenter.voxelSize/f
    ax.plot( wavelength, scatCross, label="Center", color="#4daf4a")

    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Scattering cross section, \$\sigma_s\$ (\$\SI{}{\micro\meter\squared}\$)")
    ax.legend(loc="best", frameon=False, labelspacing=0.05)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")
    plt.show()

if __name__ == "__main__":
    main()
