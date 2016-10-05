import numpy as np

class Axis:
    def __init__(self):
        self.nbins = 0
        self.min = 0.0
        self.max = 0.0

class Histogram2D:
    def __init__(self):
        self.ncol = 0
        self.nrow = 0
        self.phi = Axis()
        self.alpha = Axis()
        self.data = None
        self.deg = 180.0/np.pi
        self.unit = "rad"

    def parse(self, fname):
        readPhiData = False
        readAlphaData = False
        with open(fname, 'r') as infile:
            for line in infile:
                if ( line.find("# axis-0") != -1):
                    readPhiData = True
                    continue
                elif( line.find("# axis-1") != -1):
                    readAlphaData = True
                    continue
                elif ( line.find("# data") != -1 ):
                    break

                if ( readPhiData ):
                    self.readAxisData(line, self.phi )
                    readPhiData = False
                elif ( readAlphaData ):
                    self.readAxisData(line, self.alpha )
                    readAlphaData = False
            self.data = np.loadtxt(infile)

        self.data = self.data.reshape((self.phi.nbins, self.alpha.nbins))

    def readAxisData(self, line, axisObj):
        listContent = line.split(",")
        axisObj.nbins = int( listContent[1] )
        axisObj.min = float( listContent[2] )
        axisObj.max = float( listContent[3][:-2] )

    def display( self ):
        print ("===== Phi =====")
        print ("Nbins: %d"%(self.phi.nbins))
        print ("Min: %.2E"%(self.phi.min))
        print ("Max: %.2E"%(self.phi.max))
        print ("===== Alpha =====")
        print ("Nbins: %d"%(self.alpha.nbins))
        print ("Min: %.2E"%(self.alpha.min))
        print ("Max: %.2E"%(self.alpha.max))

    def radToDeg(self):
        if ( self.unit == "rad" ):
            self.phi.min *= self.deg
            self.phi.max *= self.deg
            self.alpha.min *= self.deg
            self.alpha.max *= self.deg
            self.unit = "deg"
