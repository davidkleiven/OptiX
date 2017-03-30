import sys
import numpy as np
from matplotlib import pyplot as plt
import subprocess
import h5py as h5

class BlochRegion:
    def __init__( self ):
        self.k = []
        self.omega = []

class BlochRun:
    def __init__(self):
        self.Ez = []
        self.temporaryFname = "temporaryEz.txt"
        self.resultFile = "results.csv"
        self.Qmax = 10       # Default harminv value
        self.modeError = 0.1 # Default harminv value
        self.fmin = 0.0
        self.fmax = 0.0
        self.amplitudeRelMax = 0.1
        self.dt = 1.0

        self.freq = []
        self.decay = []
        self.Q = []
        self.amplitude = []
        self.phase = []
        self.error = []

        self.bandDiag = plt.figure()
        self.bandAx = self.bandDiag.add_subplot(1,1,1)
        self.blochRegs = []

    def extractResonantModes( self ):
        # Export a temporary txt file with white space separated items
        assert( len(self.Ez) > 0 )
        np.savetxt( self.temporaryFname, self.Ez, delimiter=" ")

        infile = open( self.temporaryFname, 'r')
        outfile = open( self.resultFile, 'w')

        # Analyse the results with harminv
        subprocess.call(["harminv", "-t %.9f"%(self.dt), "%.9f-%.9f"%(self.fmin,self.fmax), "-Q 1000"], stdin=infile, stdout=outfile)
        outfile.close()
        infile.close()
        self.parseOutput()
        self.clean()

    def modesFFT( self ):
        values = np.zeros(10*len(self.Ez))
        values[:len(self.Ez)] = self.Ez
        freq = np.fft.fftshift( np.fft.fftfreq(len(values), d=self.dt) )
        ft = np.fft.fft(values)
        ft = np.fft.fftshift(ft)
        ft = ft[freq>self.fmin]
        freq = freq[freq>self.fmin]
        ft = ft[freq<self.fmax]
        freq = freq[freq<self.fmax]
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot( freq, np.abs(ft)**2, color="black")
        ax.set_yscale("log")

    def filterData( self, data ):
        tempFreq = data[:,0]
        data = data[tempFreq>self.fmin,:]
        tempFreq = data[:,0]
        data = data[tempFreq<self.fmax,:]

        amp = data[:,3]
        data = data[amp>self.amplitudeRelMax*amp.max(),:]
        return data

    def parseOutput( self ):
        data = np.loadtxt( self.resultFile, delimiter=",", skiprows=1)
        data = self.filterData(data)
        self.freq.extend(data[:,0])
        self.decay.extend(data[:,1])
        self.Q.extend(data[:,2])
        self.amplitude.extend(data[:,3])
        self.phase.extend(data[:,4])
        self.error.extend(data[:,5])

    def clearHarminvData( self ):
        self.freq = []
        self.decay = []
        self.Q = []
        self.amplitude = []
        self.phase = []
        self.error = []

    def clean( self ):
        subprocess.call(["rm", self.temporaryFname])
        subprocess.call(["rm", self.resultFile])

    def printInfo( self ):
        print ("Frequencies:")
        print (self.freq)
        print ("Decay rate:")
        print (self.decay)
        print ("Q-factor:")
        print (self.Q)
        print ("Amplitude:")
        print (self.amplitude)
        print ("Error:")
        print (self.error)

    def buildBlochRegs( self, hfile ):
        dset = "Run0"
        counter = 1
        while( dset in hfile.keys() ):
            self.clearHarminvData()
            group = hfile.get(dset)
            self.Ez = np.array( group.get("Ez") )
            self.fmin = group.attrs["fmin"]
            self.fmax = group.attrs["fmax"]
            self.extractResonantModes()
            breg = int(group.attrs["blochPath"])
            while (  breg >= len(self.blochBlochRegs) ):
                self.blochRegs.append(BlochRegion())
            k = np.sqrt( group.attrs["kx"]**2 + group.attrs["ky"]**2 )
            self.blochRegs[breg].k.append(k)
            self.omega.append( self.freq )
            dset = "Run%d"%(counter)
            counter += 1

    def plotBands( self ):
        assert( len(self.blochRegs) > 0 )
        klastEnd = 0.0
        colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999"]
        for reg in self.blochRegs:
            for i in range(0, len(reg.k) ):
                for j in range(0,len(reg.omega[i])):
                    self.bandAx.plot( reg.k[i], reg.omega[i][j], 'o', color=colors[j%len(colors)] )



def main( argv ):
    if ( len(argv) != 1 ):
        print ("Usage: python bandDiagram.py h5file.h5")
        return 1

    analyser = BlochRun()
    with h5.File(argv[0],'r') as hf:
        analyser.Ez = np.array( hf.get("Ez") )
        freq = np.array(hf.get("freq"))
        analyser.fmin = freq[0]/(2.0*np.pi)
        analyser.fmax = (freq[0] + freq[1]*freq[2])/(2.0*np.pi)
        analyser.dt = np.array(hf.get("dt"))[0]
    analyser.extractResonantModes()
    analyser.printInfo()
    analyser.modesFFT()
    plt.show()

if __name__ == "__main__":
    main( sys.argv[1:] )
