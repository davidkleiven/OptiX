import grazing as graz
import grazingTransmission as gzt
import numpy as np
import matplotlib
import mplLaTeX as ml
matplotlib.rcParams.update( ml.params )
from matplotlib import pyplot as plt
import colorScheme as cs
import transformDetector as td

class CombinedDetector:
    def __init__(self):
        self.n2 = 0.992+0.002j
        self.thick = 2.0
        self.reflected = graz.GrazingHandler(True)
        self.transmitted = gzt.GrazingTransmissionHandler()
        self.reflected.setEpsilonSubst(self.n2**2)
        self.transmitted.setEpsilonSubst(self.n2**2)
        self.reflected.setFilmThickness(self.thick)
        self.transmitted.setFilmThickness(self.thick)
        self.reflected.detectorTransform = td.DetectorCenterBeam()
        self.transmitted.detectorTransform = td.DetectorCenterBeam()

    def plotTotal(self, angles):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        self.reflected.prepareDWBA(angles[0])
        self.transmitted.prepareDWBA(angles[0])
     
        for i in range(0, len(angles)):
            if ( angles[i] < 1.0 ):
                label = "$\\alpha = %.1f\\alpha_c$"%(angles[i])
            elif ( angles[i] == 1 ):
                label = "$\\alpha = \\alpha_c$"%(angles[i])
            else:
                label = "$\\alpha=%d\\alpha_c$"%(angles[i])
            scatteringAngles = np.linspace(-89.0*np.pi/(180.0*self.reflected.alpha_c), 85.0*np.pi/(180.0*self.reflected.alpha_c), 100001)
            validPlus = self.reflected.adjustX(angles[i], scatteringAngles)
            scatteringAngles = np.linspace(-89.0*np.pi/(180.0*self.reflected.alpha_c), 85.0*np.pi/(180.0*self.reflected.alpha_c), 100001)
            validMinus = self.transmitted.adjustX(angles[i], scatteringAngles)

            self.reflected.prepareDWBA(angles[i])
            self.transmitted.prepareDWBA(angles[i]) 
            #alpha_f_plus = np.arctan(self.reflected.x/self.reflected.detectorPosition)/self.reflected.alpha_c
            #alpha_f_minus = np.arctan(self.transmitted.x/self.transmitted.detectorPosition)/self.transmitted.alpha_c
         
            totPluss = np.abs(self.reflected.bornTotal())**2
            firstBornPluss = np.abs(self.reflected.f1)**2
            firstBornMinus = np.abs(self.transmitted.born())**2
            totMinus = np.abs(self.transmitted.total())**2
            print len(validPlus), len(totPluss), len(firstBornPluss)
            print len(validMinus), len(totMinus), len(firstBornMinus)
        
            #scatPlus = self.reflected.detectorTransform.scatteringAngle(angles[i], alpha_f_plus)
            #scatMinus = self.reflected.detectorTransform.scatteringAngle(angles[i], alpha_f_minus)
            ax.plot( validPlus, totPluss, color=cs.COLORS[i], label=label)
            ax.plot( validMinus, totMinus, color=cs.COLORS[i])
            if ( i == len(angles)-1 ):
                ax.plot( validPlus, firstBornPluss, lw=0.3, color=cs.COLORS[len(angles)], label="BA")
            else:
                ax.plot( validPlus, firstBornPluss, lw=0.3, color=cs.COLORS[len(angles)])
            ax.plot( validMinus, firstBornMinus, lw=0.3, color=cs.COLORS[len(angles)])

        fname = "Figures/dwbaTotalPattern.pdf" 
        ax.set_xlabel(self.reflected.detectorTransform.axisLabel())
        ax.set_ylabel("Intensity (a.u.)")
        ax.set_yscale("log")
        ax.set_ylim(bottom=1E-6)
        ax.legend(loc="upper left", frameon=False, labelspacing=0.05)
        fig.savefig(fname, bbox_inches="tight")
        print ("Figure written to %s"%(fname)) 
        

# Add a main function for plotting the total spectrum
def main():
    angles = np.array([0.5,1.0,2.0,5.0])
    try:
        combHandler = CombinedDetector()
        combHandler.plotTotal(angles)
    except Exception as e:
        print str(e)

if __name__ == "__main__":
    main()
