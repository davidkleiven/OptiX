import sys
sys.path.append("../FresnelFDTD")
sys.path.append("../")
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
import colorScheme as cs
from scipy import signal

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

        freq = np.abs(freq)
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

def plot1DModes( hffile ):
    keys = hffile.keys()
    nkeys = int( len(keys)/5 )

    print ("Collecting information of frequency range")
    freqMax = -np.inf
    freqMin = np.inf
    for i in range(0, nkeys):
        freqname = "freq%d"%(i)
        freq = np.array( hffile.get(freqname) )
        freq = np.abs(freq)
        if ( np.max(freq) > freqMax ):
            freqMax = np.max(freq)
        if ( np.min(freq) < freqMin ):
            freqMin = np.min(freq)

    nFreq = 1000
    count = np.zeros(nFreq)
    amplitude = np.zeros((nkeys,nFreq))
    decay = np.zeros((nkeys,nFreq))
    bins = np.linspace(freqMin, freqMax, nFreq)
    print ("Creating histograms")
    for i in range(0,nkeys):
        freqname = "freq%d"%(i)
        decayname = "decay%d"%(i)
        ampname = "amplitude%d"%(i)
        currentAmp = np.array( hffile.get(ampname) )
        currentFreq = np.array( hffile.get(freqname) )
        currentL = 1.0/np.abs( np.array( hffile.get(decayname) ) )
        currentFreq = np.abs(currentFreq)
        localCounter = np.zeros(len(count))
        df = (freqMax-freqMin)/(nFreq-1)
        for j in range(0,len(currentFreq)):
            indx = int( (currentFreq[j]-freqMin)/df )
            count[indx] += 1
            amplitude[i,indx] = localCounter[indx]*amplitude[i,indx]+currentAmp[j]
            decay[i,indx] += localCounter[indx]*decay[i,indx]+currentL[j]
            localCounter[indx] += 1
            amplitude[i,indx] /= localCounter[indx]
            decay[i,indx] /= localCounter[indx]

    #amplitude, decay = mergeTransverse( amplitude, decay, 2)
    # Plot the 4 largest
    print ("Creating plots")
    nModes = 4
    largestIndx = np.argsort(count)[::-1][:nModes]
    largestIndx = np.sort(largestIndx)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(1,1,1)
    for i in range(0, nModes):
        x = np.linspace(0.0, 100.0, len(amplitude[:,0]))
        wLen = len(x)/10
        if ( wLen%2 == 0 ):
            wLen += 1
        inverseWav = bins[largestIndx[i]]+0.5*df
        ampSmooth = signal.savgol_filter( amplitude[:,largestIndx[i]], wLen, 3)
        ax.plot( x, ampSmooth, label="%d mm$^{-1}$"%(bins[largestIndx[i]]*1E6), color=cs.COLORS[i])

        wLen = len(x)/8
        if ( wLen%2 == 0 ):
            wLen += 1
        decaySmooth = signal.savgol_filter( decay[:,largestIndx[i]], wLen, 3)
        ax2.plot( x, decaySmooth/1E6, label="%d mm$^{-1}$"%(bins[largestIndx[i]]*1E6), color=cs.COLORS[i])
    ax.set_xlabel("Transverse position (nm)")
    ax.set_ylabel("Amplitude (a.u.)")
    ax.set_yscale("log")
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(top=10.0*ymax)
    ax.legend(loc="upper center", frameon=False, ncol=2)
    ax2.set_xlabel("Transverse position (nm)")
    ax2.set_ylabel("Decay length (mm)")
    ax2.set_yscale("log")
    ymin, ymax = ax2.get_ylim()
    ax2.set_ylim(top=10.0*ymax)
    ax2.legend(loc="upper center", frameon=False, ncol=2)
    fname = "Figures/harminvAmplitude1D.pdf"
    fig.savefig(fname, bbox_inches="tight")
    print ("Figure written to %s"%(fname))

    fname = "Figures/harminvDecay1D.pdf"
    fig2.savefig(fname, bbox_inches="tight")
    print ("Figure written to %s"%(fname))

    # Plot frequency distribution
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot((bins+0.5*df)*1E6, count, color="black", ls="steps")
    ax.set_xlabel("Inverse wavelength (mm$^{-1}$)")
    ax.set_xlim(left=-0.1*bins[-1]*1E6)
    ax.set_ylabel("Counts")
    fname = "Figures/freqDistribution.pdf"
    fig.savefig(fname, bbox_inches="tight")
    print ("Figure written to %s"%(fname))

def mergeTransverse( amplitude, decay, nMerge ):
    nrows = int( amplitude.shape[0]/nMerge )
    ncol = amplitude.shape[1]
    newAmplitude = np.zeros((nrows, ncol))
    newDecay = np.zeros((nrows,ncol))
    for i in range(0, nrows):
        newAmplitude[i,:] = np.mean(amplitude[i:i+nMerge,:], axis=0)
        newDecay[i,:] = np.mean(decay[i:i+nMerge,:], axis=0)
    return newAmplitude, newDecay

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
        #amplitudeView(hf)
        plot1DModes(hf)

if __name__ == "__main__":
    main( sys.argv[1:])
