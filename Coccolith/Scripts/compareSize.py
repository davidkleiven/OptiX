import plotSpectrum as psp
import numpy as np
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["axes.unicode_minus"] = False
mpl.rcParams["font.size"] = 18
from matplotlib import pyplot as plt

def main():
    fnameLarge = "data/CaCO3CrossSectionLong.h5"
    fnameSmall = "data/CaCO3Coccolith9ScatCross.h5"

    spLarge = psp.initSpectrum( fnameLarge )
    spSmall = psp.initSpectrum( fnameSmall )

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    scatCross = spLarge.getScatteringCrossSecion()
    f = np.linspace( spLarge.freqmin, spLarge.freqmax(), len( scatCross) )
    wavelength = spLarge.voxelSize/f
    wavemax = 800.0
    scatCross = scatCross[wavelength < wavemax]
    wavelength = wavelength[wavelength<wavemax]
    ax.plot( wavelength, scatCross, color="#e41a1c", label="A")

    scatCross = spSmall.getScatteringCrossSecion()
    f = np.linspace( spSmall.freqmin, spSmall.freqmax(), len( scatCross) )
    wavelength = spSmall.voxelSize/f
    wavemax = 800.0
    scatCross = scatCross[wavelength < wavemax]
    wavelength = wavelength[wavelength<wavemax]
    ax.plot( wavelength, scatCross, color="#377eb8", label="B")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Scattering cross section (\$\SI{}{\micro\meter\squared}\$)")
    ax.legend( loc="best", frameon=False, labelspacing=0.05)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    plt.show()

if __name__ == "__main__":
    main()
