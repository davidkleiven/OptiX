import sys
import numpy as np
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.size"] = 18
mpl.rcParams["axes.unicode_minus"] = False
from matplotlib import pyplot as plt
import extractUnitCell as euc

def main( argv ):
    if ( len(argv) != 1 ):
        print ("Usage: python extractUnitCell.py <fname_res_Nx_Ny_Nz.raw>")
        return 1
    data = np.fromfile(argv[0], dtype=np.uint8)
    res,Nx,Ny,Nz = euc.getSize(argv[0])
    data = data.reshape((Nx,Ny,Nz), order="F")
    proj = data.sum(axis=1).T
    proj = proj*255/int(proj.max())
    proj -= np.mean(proj)

    ft = np.fft.fft2( proj )
    ft = np.fft.fftshift(ft)

    kxmax = np.pi*ft.shape[0]*10/res
    kymax = np.pi*ft.shape[1]*10/res
    fig = plt.figure()
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    ax1.imshow(proj, cmap="bone", aspect="equal", extent=[0,ft.shape[0]*res/10,0,ft.shape[1]*res/10])
    maxval=np.max(np.abs(ft))
    minval=0.01*np.max(np.abs(ft))
    ax2.imshow(np.abs(ft), cmap="gist_heat", aspect="equal", norm=mpl.colors.LogNorm(minval,maxval), extent=[-kxmax,kxmax,-kymax,kymax])
    ax1.set_xlabel("$x$ (nm)")
    ax1.set_ylabel("$y$ (nm)")
    ax2.set_xlabel("$G_x (nm\$^{-1}\$)")
    ax2.set_ylabel("$G_y (nm\$^{-1}\$)")
    plt.show()

if __name__ == "__main__":
    main( sys.argv[1:] )
