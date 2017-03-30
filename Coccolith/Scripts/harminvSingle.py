import sys
import numpy as np
from matplotlib import pyplot as plt
import subprocess
import h5py as h5

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
