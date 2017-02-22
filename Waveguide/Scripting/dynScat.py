import sys
import numpy as np
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["axes.unicode_minus"] = False
mpl.rcParams["font.size"] = 18
from matplotlib import pyplot as plt
import dynScatGeo as dsg
from scipy import ndimage

def computeScatteringTerm( geometry, term ):
    assert( term >= 1 )
    x = np.linspace(geometry.xmin, geometry.xmax, geometry.Nx )
    y = np.linspace( geometry.ymin, geometry.ymax, geometry.Ny )
    X,Y = np.meshgrid(x,y)
    proj = geometry.getProjection( X, Y )
    del X, Y

    proj = np.fft.fft2( proj**term, s=(2048,2048) )
    proj = np.fft.fftshift( proj )
    return proj


def computeTotalScattering( geometry, nTerms ):
    for i in range(1, nTerms+1):
        if ( i == 1 ):
            scat = computeScatteringTerm( geometry, i )*(-1j*geometry.k)**i
        else:
            scat += (-geometry.k*1j)**i *computeScatteringTerm( geometry, i)/i
    return np.abs(scat)**2

def main():
    geo = dsg.Sphere()
    geo.k = 40.0 # nm^-1
    geo.R = 2000.0 # nm
    geo.xmin = -geo.R
    geo.xmax = geo.R
    geo.ymin = -geo.R
    geo.ymax = geo.R

    #dx = (geo.xmax-geo.xmin)/Nx
    #dy = (geo.ymax - geo.ymin)/Ny

    # Create plot illustrating the effect of the farfield as kR is increase
    dkR = [0.01, 0.1, 1.0, 10.0]
    labs=["\$10^{-2}\$", "\$10^{-1}\$", "1", "10"]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3"]
    #colors=["#1b9e77", "#d95f02", "#7570b3", "#e7298a"]
    qRmax = 10.0
    for i in range(0, len(dkR) ):
        geo.R = dkR[i]/(geo.delta*geo.k)
        geo.xmin = -geo.R
        geo.xmax = geo.R
        geo.ymin = -geo.R
        geo.ymax = geo.R
        dx = (geo.xmax-geo.xmin)/geo.Nx
        if ( i == 0 ):
            form = np.abs( geo.k*computeScatteringTerm( geo, 1 ) )**2
            com = ndimage.measurements.center_of_mass( form )
            formP = form[int(com[0]),:]
            q = np.fft.fftfreq( len(formP), d=dx)
            q = np.fft.fftshift(q)
            formP = formP[np.abs(q*geo.R) < qRmax]
            q = q[np.abs(q*geo.R)<qRmax]
            ax.plot( q*geo.R, formP, color="black", lw=4, label="\$F(q)\$")

        ff = geo.farField()
        com = ndimage.measurements.center_of_mass( ff )
        ffP = ff[int(com[0]),:]

        q = np.fft.fftfreq( len(ffP), d=dx)
        q = np.fft.fftshift(q)
        ffP = ffP[np.abs(q*geo.R) < qRmax]
        q = q[np.abs(q*geo.R)<qRmax]
        ax.plot( q*geo.R, ffP, label=labs[i], color=colors[i] )

    ax.set_yscale("log")
    ax.set_xlabel("\$qR\$")
    ax.set_ylabel("Intensity (a.u.)")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")
    ax.legend(loc="lower center", frameon=False, ncol=2 )
    plt.show()

if __name__ == "__main__":
    main()
