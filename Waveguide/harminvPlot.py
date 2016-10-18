import sys
sys.path.append("../FresnelFDTD")
import mplLaTeX as ml
import matplotlib as mpl
#mpl.rcParams.update(ml.params)
CALL_OVER_SSH = False
if ( CALL_OVER_SSH ):
    mpl.use("Agg")
import numpy as np
import h5py as h5
import json
from matplotlib import pyplot as plt
from scipy import interpolate
import transmission as trans

def amplitudeView( hffile ):
    keys = hffile.keys()
    nkeys = len(keys)/5

    fig = plt.figure()
    figdecay = plt.figure()
    ax = fig.add_subplot(1,1,1)
    axdecay = figdecay.add_subplot(1,1,1)
    cb = None
    cmin = 0.0
    cmax = 1.0
    width = 100.0
    cbDecay = None
    cminDec = 1E-2
    cmaxDec = 1E3
    for i in range(0, nkeys):
        amp = np.array( hffile.get("amplitude%d"%(i)) )
        freq = np.array( hffile.get("freq%d"%(i)))
        decay = np.array( hffile.get("decay%d"%(i)))
        x = np.zeros(len(amp))+width*i/nkeys
        # Filter out negative frequencies
        #x = x[freq>0.0]
        #amp = amp[freq>0.0]
        #freq = freq[freq>0.0]

        # Filter out very low amplitudes
        #x = x[np.abs(amp)>0.05*np.max(amp)]
        #freq = freq[np.abs(amp)>0.05*np.max(amp)]
        #amp = amp[np.abs(amp)>0.05*np.max(amp)]

        freq *= 1E6
        L = np.abs(1.0/decay)/1E6
        s = ax.scatter(freq, x, c=np.abs(amp), cmap="coolwarm", s=0.1, lw=0 )
        sdec = axdecay.scatter(freq, x, c=L, cmap="coolwarm", s=0.1, lw=0, norm=mpl.colors.LogNorm(vmin=cminDec,vmax=cmaxDec))
        if ( cb is None ):
            cb = fig.colorbar(s)

        if ( cbDecay is None ):
            cbDecay = figdecay.colorbar(sdec)

        if ( np.max(L) > cmaxDec ):
            cmaxDec = np.max(L)
        if ( np.min(L) < cminDec ):
            cminDec = np.min(L)
        s.set_clim([cmin,cmax])
        #sdec.set_clim([cminDec, cmaxDec])

    print ("Longest decay length: %.2E mm"%(cmaxDec))
    ax.set_xlabel("Inverse wavelength (mm$^{-1}$)")
    axdecay.set_xlabel("Inverse wavelength (mm$^{-1}$)")
    ax.set_ylabel("Transverse position (nm)")
    axdecay.set_ylabel("Transverse position (nm)")
    xmin = -1E-4*1E6
    xmax = 1E-4*1E6
    ax.set_xlim(xmin, xmax)
    axdecay.set_xlim(xmin,xmax)
    #ax.set_xscale("log")
    fname = "Figures/harminv.jpeg"
    fig.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))

    fname = "Figures/harminv_decay.jpeg"
    figdecay.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))

def amplitudeFFT( hffile ):
    data = np.array( hffile.get( "dataset" ))

    extent = [0.0, 100.0, 0.0, 2.0*np.pi/0.5]
    plt.imshow( np.abs(data)**2, extent=extent, aspect=1.0, cmap="coolwarm")#, norm=mpl.colors.LogNorm())
    plt.xlabel("$k_z$")
    plt.ylabel("$x$ (nm)")
    plt.gca().set_aspect((extent[1]-extent[0])/(extent[3]-extent[2]))
    plt.colorbar()
    fname = "Figures/fieldFFT.jpeg"
    plt.savefig(fname, bbox_inches="tight", dpi=800)
    print ("Figure written to %s"%(fname))

def plot1D( hffile ):
    print hffile.keys()
    mpl.rcParams.update(mpl.rcParamsDefault)
    data = np.array( hffile.get("transformedData0") )
    plt.plot(np.log(np.abs(data)), 'k')
    plt.show()

def main(argv):
    for arg in argv:
        if ( arg.find("--file=") != -1 ):
            fname = arg.split("--file=")[1]
        else:
            print ("Unknown argument %s"%(arg))
            return 0

    with h5.File(fname, 'r') as hf:
        #plot1D(hf)
        #amplitudeFFT(hf)
        amplitudeView(hf)

if __name__ == "__main__":
    main( sys.argv[1:])
